// ConsoleApplication.cpp : This file contains the 'main' function. Program
// execution begins and ends there.
//
#include <imstkParallelUtils.h>
#include <imstkParallelFor.h>

#include <string.h>
#include <chrono>
#include <iostream>
#include <string>
#include <thread>

#include "HardwareAPI.h"
#include "imstkMath.h"

#include "imstkHaplyDeviceClient.h"
#include "imstkHaplyDeviceManager.h"

#include "imstkOpenHapticDeviceClient.h"
#include "imstkOpenHapticDeviceManager.h"

#include "imstkSimulationManager.h"
#include "imstkSceneManager.h"
#include "imstkScene.h"

#include <iostream>

#include <HardwareAPI.h>

using namespace imstk;

class MyHandle : public Haply::HardwareAPI::Devices::Handle
{
public:
    MyHandle(std::iostream* stream) : Handle(stream)
    {
        calibrated       = false;
        quaternionOffset = new float[4]{ 1, 0, 0, 0 };
        statusResponse   = Haply::HardwareAPI::Devices::Handle::HandleStatusResponse();
        statusResponseCalibrated = Haply::HardwareAPI::Devices::Handle::HandleStatusResponse();
    }

    //  Devices::Handle::HandleStatusResponse* debug() {
    //return &statusResponse;
    //  };

    void Calibrate()
    {
        float knownQuaternion[4] = { 1, 0, 0, 0 };
        // Copy the current quaternion
        float currentQuaternion[4] = { statusResponse.quaternion[0], statusResponse.quaternion[1], statusResponse.quaternion[2], statusResponse.quaternion[3] };

        // Determine the transformation required to make the current quaternion equal to the known Quaternion
        QuaternionConjugate(currentQuaternion);
        QuaternionMultiply(knownQuaternion, currentQuaternion);

        // Store the quaternion offset
        quaternionOffset[0] = knownQuaternion[0];
        quaternionOffset[1] = knownQuaternion[1];
        quaternionOffset[2] = knownQuaternion[2];
        quaternionOffset[3] = knownQuaternion[3];

        calibrated = true;
    }

    bool IsCalibrated() { return calibrated; }
    Haply::HardwareAPI::Devices::Handle::HandleInfoResponse GetInfoResponse() { return infoResponse; }
    Haply::HardwareAPI::Devices::Handle::HandleStatusResponse GetStatusResponse()
    {
        if (calibrated)
        {
            statusResponseCalibrated = statusResponse;

            // Create quaternion offset copy to modify the value in the function
            float quaternionOffsetCopy[4] = { quaternionOffset[0], quaternionOffset[1], quaternionOffset[2], quaternionOffset[3] };
            float quaternion[4] = { statusResponse.quaternion[0], statusResponse.quaternion[1], statusResponse.quaternion[2], statusResponse.quaternion[3] };

            // Apply Quaternion calibration transformation to the current quaternion
            QuaternionMultiply(quaternionOffsetCopy, quaternion);

            // Assign the calibrated quaternion to the status response
            statusResponseCalibrated.quaternion[0] = quaternionOffsetCopy[0];
            statusResponseCalibrated.quaternion[1] = quaternionOffsetCopy[1];
            statusResponseCalibrated.quaternion[2] = quaternionOffsetCopy[2];
            statusResponseCalibrated.quaternion[3] = quaternionOffsetCopy[3];

            // Return status Response
            return statusResponseCalibrated;
        }
        else
        {
            return statusResponse;
        }
    }

    Haply::HardwareAPI::Devices::Handle::HandleErrorResponse GetErrorResponse() { return errorResponse; }

protected:
    void OnReceiveHandleInfo(HandleInfoResponse& response) override { infoResponse = response; }
    void OnReceiveHandleStatusMessage(HandleStatusResponse& response) override { statusResponse = response; }
    void OnReceiveHandleErrorResponse(HandleErrorResponse& response) override { errorResponse = response; }

private:
    bool calibrated;

    Haply::HardwareAPI::Devices::Handle::HandleInfoResponse   infoResponse;
    Haply::HardwareAPI::Devices::Handle::HandleStatusResponse statusResponse;
    Haply::HardwareAPI::Devices::Handle::HandleStatusResponse statusResponseCalibrated;
    Haply::HardwareAPI::Devices::Handle::HandleErrorResponse  errorResponse;

    void QuaternionMultiply(float* q1, float* q2)
    {
        float q[3];
        q[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
        q[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
        q[2] = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
        q[3] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];

        q1[0] = q[0];
        q1[1] = q[1];
        q1[2] = q[2];
        q1[3] = q[3];
    }

    void QuaternionConjugate(float* q)
    {
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
    }

    float* quaternionOffset;
};

int
test1()
{
    Logger::startLogger();

    char* portName;

    std::string portNames[256];
    auto        nbport =
        Haply::HardwareAPI::Devices::DeviceDetection::AutoDetectInverse3(
            portNames);
    printf("Found %d ports\n", nbport);
    if (nbport > 0)
    {
        int index = nbport - 1;
        portName = _strdup(portNames[index].c_str());
    }
    else
    {
        std::cout << "No Inverse3 found" << std::endl;
        return -1;
    }

    printf("Using port %s\n", portName);
    Haply::HardwareAPI::IO::SerialStream  serial_stream(portName, false);
    char                                  open_return = serial_stream.OpenDevice();
    Haply::HardwareAPI::Devices::Inverse3 inverse3(&serial_stream);

    MyHandle*                             handle;
    Haply::HardwareAPI::IO::SerialStream* handleStream;
    bool                                  handleFound = false;
    // Find Handle
    std::string handlePorts[100];
    int         handleCount = Haply::HardwareAPI::Devices::DeviceDetection::AutoDetectHandle(handlePorts);
    // If handle is found, create a new handle object
    if (handleCount == 1)
    {
        printf("Handle Found");
        char* portName = _strdup(handlePorts[0].c_str());

        // Create a new stream and a new device
        handleStream = new Haply::HardwareAPI::IO::SerialStream(portName);
        handle       = new MyHandle(handleStream);
        handleFound  = true;
    }
    // Initialize variables used inside of the loop
    Haply::HardwareAPI::Devices::Inverse3::EndEffectorStateResponse ee_response;
    Haply::HardwareAPI::Devices::Handle::HandleStatusResponse       handle_status;
    unsigned char                                                   lastReturnType;

    // Wake up the handle. If the handle is awake, the first response may not be correct response type, so we wait until
    if (handleFound)
    {
        handle->SendDeviceWakeup();
        while (lastReturnType != 0xD0)
        {
            handle->Receive(&lastReturnType);
        }
        Haply::HardwareAPI::Devices::Handle::HandleInfoResponse handleInfo = handle->GetInfoResponse();
        // Calibrate the handle. Make sure to point the handle toward the screen
        handle->Calibrate();
    }

    inverse3.SendDeviceWakeup();
    inverse3.ReceiveDeviceInfo();

    /* std::cout << std::endl << "Press ENTER to continue . . .";
     std::cin.get();*/

    Vec3f pos   = Vec3f::Zero();
    Vec3f vel   = Vec3f::Zero();
    Vec3f force = Vec3f::Zero();

    // Bounds ~ -0.5, 0.5
    while (true)
    {
        /*ee_response = inverse3.EndEffectorForce({ force[0], force[1], force[2] });
        printf("pos: %f, %f, %f\n", ee_response.position[0], ee_response.position[1], ee_response.position[2]);*/

        inverse3.SendEndEffectorForce(force.data());
        inverse3.ReceiveEndEffectorState(pos.data(), vel.data());

        // Get handle update
        if (handleFound)
        {
            handle->Receive(&lastReturnType);
            // Check the handle message type and print handle states
            if (lastReturnType == HAPLY_TOOL_STATUS)
            {
                handle_status = handle->GetStatusResponse();
            }
        }

        double planePos = 0.1;
        double ks       = 250.0;
        if (pos[2] < planePos)
        {
            force[2] = (planePos - pos[2]) * ks;
        }
        else
        {
            force[2] = 0.0;
        }
        //printf("Applying force: %f\n", force[2]);
    }
}

void
test2()
{
    auto manager = std::make_shared<HaplyDeviceManager>();
    auto client  = manager->makeDeviceClient();

    manager->init();

    while (true)
    {
        manager->update();

        const Vec3d& pos      = client->getPosition();
        double       planePos = 0.1;
        double       ks       = 500.0;
        if (pos[1] < planePos)
        {
            client->setForce(Vec3d(0.0, (planePos - pos[1]) * ks, 0.0));
        }
        else
        {
            client->setForce(Vec3d::Zero());
        }
    }

    manager->uninit();
}

int
main(int argc, char* argv[])
{
    Logger::startLogger();

    test2();

    //test1();
}