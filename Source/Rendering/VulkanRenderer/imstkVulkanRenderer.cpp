/*=========================================================================

   Library: iMSTK

   Copyright (c) Kitware, Inc. & Center for Modeling, Simulation,
   & Imaging in Medicine, Rensselaer Polytechnic Institute.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

=========================================================================*/

#include "imstkVulkanRenderer.h"

namespace imstk
{
VulkanRenderer::VulkanRenderer(std::shared_ptr<Scene> scene)
{
    m_scene = scene;
}

void
VulkanRenderer::initialize(unsigned int width, unsigned int height)
{
    m_width = width;
    m_height = height;
    // If debug mode, enable validation layer (slower performance)
#ifndef NDEBUG
    m_layers.push_back(VulkanValidation::getValidationLayer());
    m_extensions.push_back(VulkanValidation::getValidationExtension());
#endif

    // Instance of a Vulkan application
    VkInstanceCreateInfo m_creationInfo;
    m_creationInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    m_creationInfo.pNext = nullptr;
    m_creationInfo.flags = 0;
    m_creationInfo.pApplicationInfo = nullptr;
    m_creationInfo.enabledLayerCount = (uint32_t)m_layers.size();
    m_creationInfo.ppEnabledLayerNames = &m_layers[0];
    m_creationInfo.enabledExtensionCount = (uint32_t)m_extensions.size();
    m_creationInfo.ppEnabledExtensionNames = &m_extensions[0];

    std::cout << "\n" << "Vulkan Renderer Information:" << std::endl;

    for (int i = 0; i < m_extensions.size(); i++)
    {
        std::cout << "Enabled extension: " << m_creationInfo.ppEnabledExtensionNames[i] << std::endl;
    }

    m_instance = new VkInstance();

    vkCreateInstance(&m_creationInfo, nullptr, m_instance);

#ifndef NDEBUG
    VkDebugReportCallbackCreateInfoEXT debugReportInfo;
    debugReportInfo.sType = VK_STRUCTURE_TYPE_DEBUG_REPORT_CALLBACK_CREATE_INFO_EXT;
    debugReportInfo.pNext = nullptr;
    debugReportInfo.flags =
        VK_DEBUG_REPORT_WARNING_BIT_EXT |
        VK_DEBUG_REPORT_ERROR_BIT_EXT;
    debugReportInfo.pfnCallback = &(VulkanValidation::debugReportCallback);
    debugReportInfo.pUserData = nullptr;

    PFN_vkCreateDebugReportCallbackEXT createCallback =
        (PFN_vkCreateDebugReportCallbackEXT)vkGetInstanceProcAddr(*m_instance, "vkCreateDebugReportCallbackEXT");
    createCallback(*m_instance, &debugReportInfo, nullptr, &m_debugReportCallback);
#endif

    auto camera = m_scene->getCamera();
    m_fov = (float)glm::radians(camera->getViewAngle());

    // Setup logical devices
    this->setupGPUs();
    this->printGPUs();

    // Setup command pool(s) - right now we just have one
    this->setupCommandPools();
    this->buildCommandBuffer();
    this->setupRenderPasses();
    this->setupSynchronization();
    this->setupMemoryManager();
    this->createGlobalUniformBuffers();
    this->createShadowMaps(m_shadowMapResolution);

    std::vector<VkPipeline> graphicsPipelines;
    std::vector<VkGraphicsPipelineCreateInfo> graphicsPipelinesInfo;

    VkPhysicalDeviceProperties deviceProperties;
    vkGetPhysicalDeviceProperties(m_renderPhysicalDevice, &deviceProperties);
    m_deviceLimits = deviceProperties.limits;

    VkPipelineCacheCreateInfo pipelineCacheCreateInfo;
    pipelineCacheCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_CACHE_CREATE_INFO;
    pipelineCacheCreateInfo.pNext = nullptr;
    pipelineCacheCreateInfo.flags = 0;
    pipelineCacheCreateInfo.initialDataSize = 0;
    pipelineCacheCreateInfo.pInitialData = nullptr;

    vkCreatePipelineCache(m_renderDevice, &pipelineCacheCreateInfo, nullptr, &m_pipelineCache);
}

void
VulkanRenderer::setupGPUs()
{
    // Prevent devices from being set up multiple times
    if (m_physicalDeviceCount != 0)
    {
        return;
    }

    // Setup physical devices
    vkEnumeratePhysicalDevices(*m_instance, &m_physicalDeviceCount, nullptr);
    m_physicalDevices = new VkPhysicalDevice[(int)m_physicalDeviceCount]();
    vkEnumeratePhysicalDevices(*m_instance, &m_physicalDeviceCount, m_physicalDevices);
    m_renderPhysicalDevice = m_physicalDevices[0];

    // Get render queue family
    vkGetPhysicalDeviceQueueFamilyProperties(m_renderPhysicalDevice, &m_queueFamilyPropertiesCount, nullptr);
    m_queueFamilyProperties = new VkQueueFamilyProperties[(int)m_queueFamilyPropertiesCount]();
    vkGetPhysicalDeviceQueueFamilyProperties(m_renderPhysicalDevice, &m_queueFamilyPropertiesCount, m_queueFamilyProperties);

    m_renderQueueFamily = 0;

    for (int i = 0; i < (int)(m_queueFamilyPropertiesCount); i++)
    {
        if (m_queueFamilyProperties[i].queueFlags & VK_QUEUE_GRAPHICS_BIT)
        {
            m_renderQueueFamily = i;
            break;
        }
    }

    std::vector<float> priorities(m_queueFamilyProperties[m_renderQueueFamily].queueCount);
    std::fill(priorities.begin(), priorities.end(), 0.0f);
    priorities[m_renderQueueFamily] = 1.0;

    //Setup logical devices
    m_deviceCount = m_physicalDeviceCount;
    m_devices = new VkDevice[(int)(m_deviceCount)]();

    VkDeviceQueueCreateInfo queueInfo;
    queueInfo.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
    queueInfo.pNext = nullptr;
    queueInfo.flags = 0;
    queueInfo.queueFamilyIndex = m_renderQueueFamily;
    queueInfo.queueCount = 1;// m_queueFamilyProperties[m_renderQueueFamily].queueCount;
    queueInfo.pQueuePriorities = &priorities[0];

    // The display system isn't part of the Vulkan core
    char * deviceExtensions[1];
    deviceExtensions[0] = "VK_KHR_swapchain";

    // Enabling optional Vulkan features
    VkPhysicalDeviceFeatures deviceFeatures;
    vkGetPhysicalDeviceFeatures(m_physicalDevices[0], &deviceFeatures);

    VkPhysicalDeviceFeatures features = {VK_FALSE};
    features.fillModeNonSolid = VK_TRUE;
    features.tessellationShader = VK_TRUE;
    features.wideLines = deviceFeatures.wideLines;

    if (features.wideLines == VK_TRUE)
    {
        m_supportsWideLines = true;
    }

    VkDeviceCreateInfo deviceInfo;
    deviceInfo.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
    deviceInfo.pNext = nullptr;
    deviceInfo.flags = 0;
    deviceInfo.queueCreateInfoCount = 1;
    deviceInfo.pQueueCreateInfos = &queueInfo;
#ifndef NDEBUG
    deviceInfo.enabledLayerCount = (uint32_t)m_layers.size();
    deviceInfo.ppEnabledLayerNames = &m_layers[0];
#else
    deviceInfo.enabledLayerCount = 0;
    deviceInfo.ppEnabledLayerNames = nullptr;
#endif
    deviceInfo.enabledExtensionCount = 1;
    deviceInfo.ppEnabledExtensionNames = deviceExtensions;
    deviceInfo.pEnabledFeatures = &features;

    for (int i = 0; i < (int)(m_physicalDeviceCount); i++)
    {
        vkCreateDevice(m_physicalDevices[i], &deviceInfo, nullptr, &m_devices[i]);
    }

    // This decision needs some work, may pick weaker device
    m_renderDevice = m_devices[0];

    // Get the first render-capable queue
    vkGetDeviceQueue(m_renderDevice, m_renderQueueFamily, 0, &m_renderQueue);
}

void
VulkanRenderer::printGPUs()
{
    this->setupGPUs();

    VkPhysicalDeviceProperties properties;

    std::cout << "Devices:" << std::endl;

    for (int i = 0; i < (int)(m_physicalDeviceCount); i++)
    {
        vkGetPhysicalDeviceProperties(m_physicalDevices[i], &properties);
        properties.limits.maxMemoryAllocationCount;
        std::cout << (i + 1) << ". " << properties.deviceName << std::endl;
    }
}

void
VulkanRenderer::setupCommandPools()
{
    // Create command pools (only one for now)
    VkCommandPoolCreateInfo commandPoolInfo;
    commandPoolInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    commandPoolInfo.pNext = nullptr;
    commandPoolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    commandPoolInfo.queueFamilyIndex = 0;

    vkCreateCommandPool(m_renderDevice, &commandPoolInfo, nullptr, &m_renderCommandPool);

    // Create command pools (only one for now)
    commandPoolInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    commandPoolInfo.pNext = nullptr;
    commandPoolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    commandPoolInfo.queueFamilyIndex = 0;

    vkCreateCommandPool(m_renderDevice, &commandPoolInfo, nullptr, &m_postProcessingCommandPool);
}

void
VulkanRenderer::buildCommandBuffer()
{
    // Build command buffer
    VkCommandBufferAllocateInfo commandBufferInfo;
    commandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    commandBufferInfo.pNext = nullptr;
    commandBufferInfo.commandPool = m_renderCommandPool;
    commandBufferInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    commandBufferInfo.commandBufferCount = m_buffering;

    m_renderCommandBuffer.resize(m_buffering);
    vkAllocateCommandBuffers(m_renderDevice, &commandBufferInfo, &m_renderCommandBuffer[0]);

    commandBufferInfo.commandPool = m_postProcessingCommandPool;
    m_postProcessingCommandBuffer.resize(m_buffering);
    vkAllocateCommandBuffers(m_renderDevice, &commandBufferInfo, &m_postProcessingCommandBuffer[0]);
}

void
VulkanRenderer::setupRenderPasses()
{
    // Number of geometry passes
    VulkanRenderPassGenerator::generateOpaqueRenderPass(m_renderDevice, m_opaqueRenderPass, m_samples);
    VulkanRenderPassGenerator::generateDecalRenderPass(m_renderDevice, m_decalRenderPass, m_samples);
    VulkanRenderPassGenerator::generateDepthRenderPass(m_renderDevice, m_depthRenderPass, m_samples);
}

void
VulkanRenderer::resizeFramebuffers(VkSwapchainKHR * swapchain, int width, int height)
{
    m_width = width;
    m_height = height;

    this->deleteFramebuffers();

    this->initializeFramebuffers(swapchain);

    //std::vector<VkPipeline> pipelines;
    //std::vector<VkGraphicsPipelineCreateInfo> pipelineInfos;

    //for (int i = 0; i < m_renderDelegates.size(); i++)
    //{
    //auto material = m_renderDelegates[i]->m_material;
    //vkDestroyPipeline(m_renderDevice, material->m_pipeline, nullptr);
    //material->initialize(this);
    //}
}

void
VulkanRenderer::initializeFramebufferImages(VkSwapchainKHR * swapchain)
{
    m_mipLevels = std::log2(std::max(m_width, m_height)) + 1;

    m_swapchain = swapchain;
    vkGetSwapchainImagesKHR(m_renderDevice, *m_swapchain, &m_swapchainImageCount, nullptr);

    // Depth image
    VkImageCreateInfo depthImageInfo;
    depthImageInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
    depthImageInfo.pNext = nullptr;
    depthImageInfo.flags = 0;
    depthImageInfo.imageType = VK_IMAGE_TYPE_2D;
    depthImageInfo.format = VK_FORMAT_D32_SFLOAT;
    depthImageInfo.extent = { m_width, m_height, 1 };
    depthImageInfo.mipLevels = 1;
    depthImageInfo.arrayLayers = 1;
    depthImageInfo.samples = m_samples;
    depthImageInfo.tiling = VK_IMAGE_TILING_OPTIMAL;
    depthImageInfo.usage = VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT
                           | VK_IMAGE_USAGE_SAMPLED_BIT;
    depthImageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    depthImageInfo.queueFamilyIndexCount = 0;
    depthImageInfo.pQueueFamilyIndices = nullptr;
    depthImageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

    m_depthImage.resize(m_mipLevels);
    m_depthImage[0] = m_memoryManager.requestImage(m_renderDevice,
                depthImageInfo,
                VulkanMemoryType::FRAMEBUFFER);

    for (uint32_t i = 1; i < m_mipLevels; i++)
    {
        depthImageInfo.usage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT;
        depthImageInfo.format = VK_FORMAT_R32_SFLOAT;
        depthImageInfo.extent = { std::max(m_width >> i, 1u), std::max(m_height >> i, 1u), 1 };
        m_depthImage[i] = m_memoryManager.requestImage(m_renderDevice,
                        depthImageInfo,
                        VulkanMemoryType::FRAMEBUFFER);
    }

    // Normal image
    auto normalImageInfo = depthImageInfo;
    normalImageInfo.extent = { m_width, m_height, 1 };
    normalImageInfo.usage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT |
                            VK_IMAGE_USAGE_SAMPLED_BIT;
    normalImageInfo.format = VK_FORMAT_R8G8B8A8_SNORM;

    m_normalImage = m_memoryManager.requestImage(m_renderDevice,
                normalImageInfo,
                VulkanMemoryType::FRAMEBUFFER);

    // HDR image
    auto HDRImageInfo = depthImageInfo;
    HDRImageInfo.format = VK_FORMAT_R16G16B16A16_SFLOAT;
    HDRImageInfo.mipLevels = 1;
    HDRImageInfo.usage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT
                         | VK_IMAGE_USAGE_SAMPLED_BIT;

    m_HDRImage[0].resize(m_mipLevels);
    m_HDRImage[1].resize(m_mipLevels);
    m_HDRImage[2].resize(m_mipLevels);
    for (uint32_t i = 0; i < m_mipLevels; i++)
    {
        HDRImageInfo.extent = { std::max(m_width >> i, 1u), std::max(m_height >> i, 1u), 1 };
        m_HDRImage[0][i] = m_memoryManager.requestImage(m_renderDevice, HDRImageInfo, VulkanMemoryType::FRAMEBUFFER);
        m_HDRImage[1][i] = m_memoryManager.requestImage(m_renderDevice, HDRImageInfo, VulkanMemoryType::FRAMEBUFFER);
        m_HDRImage[2][i] = m_memoryManager.requestImage(m_renderDevice, HDRImageInfo, VulkanMemoryType::FRAMEBUFFER);
    }

    for (uint32_t i = 0; i < m_swapchainImageCount; i++)
    {
        m_HDRTonemaps.resize(m_swapchainImageCount);
    }

    // AO image
    auto AOImageInfo = depthImageInfo;
    AOImageInfo.format = VK_FORMAT_R8_UNORM;
    AOImageInfo.mipLevels = 1;
    AOImageInfo.usage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT
                        | VK_IMAGE_USAGE_SAMPLED_BIT;
    AOImageInfo.extent = { m_width / 2, m_height / 2, 1 };

    m_halfAOImage[0] = m_memoryManager.requestImage(m_renderDevice, AOImageInfo, VulkanMemoryType::FRAMEBUFFER);
    m_halfAOImage[1] = m_memoryManager.requestImage(m_renderDevice, AOImageInfo, VulkanMemoryType::FRAMEBUFFER);

    // Create image views
    m_depthImageView.resize(m_mipLevels);

    for (uint32_t i = 0; i < m_mipLevels; i++)
    {
        VkImageSubresourceRange subresourceRange;
        subresourceRange.aspectMask = i == 0 ? VK_IMAGE_ASPECT_DEPTH_BIT : VK_IMAGE_ASPECT_COLOR_BIT;
        subresourceRange.baseMipLevel = 0;
        subresourceRange.levelCount = 1;
        subresourceRange.baseArrayLayer = 0;
        subresourceRange.layerCount = 1;

        VkImageViewCreateInfo imageViewInfo;
        imageViewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
        imageViewInfo.pNext = nullptr;
        imageViewInfo.flags = 0;
        imageViewInfo.image = *m_depthImage[i]->getImage();
        imageViewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;

        imageViewInfo.format = i == 0 ? VK_FORMAT_D32_SFLOAT : VK_FORMAT_R32_SFLOAT;

        imageViewInfo.components = VulkanDefaults::getDefaultComponentMapping();
        imageViewInfo.subresourceRange = subresourceRange;

        vkCreateImageView(m_renderDevice, &imageViewInfo, nullptr, &m_depthImageView[i]);
    }

    {
        VkImageSubresourceRange subresourceRange;
        subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        subresourceRange.baseMipLevel = 0;
        subresourceRange.levelCount = 1;
        subresourceRange.baseArrayLayer = 0;
        subresourceRange.layerCount = 1;

        VkImageViewCreateInfo imageViewInfo;
        imageViewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
        imageViewInfo.pNext = nullptr;
        imageViewInfo.flags = 0;
        imageViewInfo.image = *m_normalImage->getImage();
        imageViewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
        imageViewInfo.format = VK_FORMAT_R8G8B8A8_SNORM;
        imageViewInfo.components = VulkanDefaults::getDefaultComponentMapping();
        imageViewInfo.subresourceRange = subresourceRange;

        vkCreateImageView(m_renderDevice, &imageViewInfo, nullptr, &m_normalImageView);
    }

    {
        VkImageSubresourceRange subresourceRange;
        subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        subresourceRange.baseMipLevel = 0;
        subresourceRange.levelCount = 1;
        subresourceRange.baseArrayLayer = 0;
        subresourceRange.layerCount = 1;

        VkImageViewCreateInfo imageViewInfo;
        imageViewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
        imageViewInfo.pNext = nullptr;
        imageViewInfo.flags = 0;
        imageViewInfo.image = *m_halfAOImage[0]->getImage();
        imageViewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
        imageViewInfo.format = VK_FORMAT_R8_UNORM;
        imageViewInfo.components = VulkanDefaults::getDefaultComponentMapping();
        imageViewInfo.subresourceRange = subresourceRange;

        vkCreateImageView(m_renderDevice, &imageViewInfo, nullptr, &m_halfAOImageView[0]);

        imageViewInfo.image = *m_halfAOImage[1]->getImage();
        vkCreateImageView(m_renderDevice, &imageViewInfo, nullptr, &m_halfAOImageView[1]);
    }

    VkSamplerCreateInfo samplerInfo;
    samplerInfo.sType = VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO;
    samplerInfo.pNext = nullptr;
    samplerInfo.flags = 0;
    samplerInfo.magFilter = VK_FILTER_LINEAR;
    samplerInfo.minFilter = VK_FILTER_LINEAR;
    samplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_LINEAR; // Trilinear interpolation
    samplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_BORDER;
    samplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_BORDER;
    samplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_BORDER;
    samplerInfo.mipLodBias = 0.0;
    samplerInfo.anisotropyEnable = VK_FALSE;
    samplerInfo.maxAnisotropy = 1.0;
    samplerInfo.compareEnable = VK_FALSE;
    samplerInfo.compareOp = VK_COMPARE_OP_ALWAYS;
    samplerInfo.minLod = 0;
    samplerInfo.maxLod = 0;
    samplerInfo.borderColor = VK_BORDER_COLOR_FLOAT_TRANSPARENT_BLACK;
    samplerInfo.unnormalizedCoordinates = VK_FALSE;

    vkCreateSampler(m_renderDevice, &samplerInfo, nullptr, &m_HDRImageSampler);

    for (int i = 0; i < 3; i++)
    {
        m_HDRImageView[i].resize(m_mipLevels);

        for (uint32_t j = 0; j < m_mipLevels; j++)
        {
            VkImageSubresourceRange subresourceRange;
            subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
            subresourceRange.baseMipLevel = 0;
            subresourceRange.levelCount = 1;
            subresourceRange.baseArrayLayer = 0;
            subresourceRange.layerCount = 1;

            VkImageViewCreateInfo imageViewInfo;
            imageViewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
            imageViewInfo.pNext = nullptr;
            imageViewInfo.flags = 0;
            imageViewInfo.image = *m_HDRImage[i][j]->getImage();
            imageViewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
            imageViewInfo.format = VK_FORMAT_R16G16B16A16_SFLOAT;
            imageViewInfo.components = VulkanDefaults::getDefaultComponentMapping();
            imageViewInfo.subresourceRange = subresourceRange;

            vkCreateImageView(m_renderDevice, &imageViewInfo, nullptr, &m_HDRImageView[i][j]);
        }
    }
}

void
VulkanRenderer::initializeFramebuffers(VkSwapchainKHR * swapchain)
{
    // Get images from surface (color images)
    m_swapchain = swapchain;
    vkGetSwapchainImagesKHR(m_renderDevice, *m_swapchain, &m_swapchainImageCount, nullptr);
    m_swapchainImages.resize(m_swapchainImageCount);
    vkGetSwapchainImagesKHR(m_renderDevice, *m_swapchain, &m_swapchainImageCount, &m_swapchainImages[0]);
    m_swapchainImageViews.resize(m_swapchainImageCount);

    m_swapchainImageSamplers.resize(m_swapchainImageCount);
    for (uint32_t i = 0; i < m_swapchainImageCount; i++)
    {
        VkImageSubresourceRange subresourceRange;
        subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        subresourceRange.baseMipLevel = 0;
        subresourceRange.levelCount = 1;
        subresourceRange.baseArrayLayer = 0;
        subresourceRange.layerCount = 1;

        VkImageViewCreateInfo imageViewInfo;
        imageViewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
        imageViewInfo.pNext = nullptr;
        imageViewInfo.flags = 0;
        imageViewInfo.image = m_swapchainImages[i];
        imageViewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
        imageViewInfo.format = VK_FORMAT_B8G8R8A8_SRGB;
        imageViewInfo.components = VulkanDefaults::getDefaultComponentMapping();
        imageViewInfo.subresourceRange = subresourceRange;

        vkCreateImageView(m_renderDevice, &imageViewInfo, nullptr, &m_swapchainImageViews[i]);

        VkSamplerCreateInfo samplerInfo;
        samplerInfo.sType = VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO;
        samplerInfo.pNext = nullptr;
        samplerInfo.flags = 0;
        samplerInfo.magFilter = VK_FILTER_LINEAR;
        samplerInfo.minFilter = VK_FILTER_LINEAR;
        samplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_LINEAR; // Trilinear interpolation
        samplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_REPEAT;
        samplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_REPEAT;
        samplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_REPEAT;
        samplerInfo.mipLodBias = 0.0;
        samplerInfo.anisotropyEnable = VK_FALSE; // TODO:: add option to enable
        samplerInfo.maxAnisotropy = 1.0;
        samplerInfo.compareEnable = VK_FALSE;
        samplerInfo.compareOp = VK_COMPARE_OP_ALWAYS;
        samplerInfo.minLod = 0;
        samplerInfo.maxLod = 0;
        samplerInfo.borderColor = VK_BORDER_COLOR_FLOAT_TRANSPARENT_BLACK;
        samplerInfo.unnormalizedCoordinates = VK_FALSE;

        vkCreateSampler(m_renderDevice, &samplerInfo, nullptr, &m_swapchainImageSamplers[i]);
    }

    this->initializePostProcesses();

    m_opaqueFramebuffer =
        std::make_shared<VulkanFramebuffer>(m_memoryManager, m_width, m_height, false, m_samples);
    m_opaqueFramebuffer->setColor(&m_HDRImageView[0][0], VK_FORMAT_R16G16B16A16_SFLOAT);
    m_opaqueFramebuffer->setSpecular(&m_HDRImageView[1][0], VK_FORMAT_R16G16B16A16_SFLOAT);
    m_opaqueFramebuffer->setDepth(&m_depthImageView[0], VK_FORMAT_D32_SFLOAT);
    m_opaqueFramebuffer->setNormal(&m_normalImageView, VK_FORMAT_R8G8B8A8_SNORM);
    m_opaqueFramebuffer->initializeFramebuffer(&m_opaqueRenderPass);

    m_decalFramebuffer =
        std::make_shared<VulkanFramebuffer>(m_memoryManager, m_width, m_height, false, m_samples);
    m_decalFramebuffer->setColor(&m_HDRImageView[0][0], VK_FORMAT_R16G16B16A16_SFLOAT);
    m_decalFramebuffer->setSpecular(&m_HDRImageView[1][0], VK_FORMAT_R16G16B16A16_SFLOAT);
    m_decalFramebuffer->setDepth(&m_depthImageView[0], VK_FORMAT_D32_SFLOAT);
    m_decalFramebuffer->initializeFramebuffer(&m_decalRenderPass);

    m_depthFramebuffer =
        std::make_shared<VulkanFramebuffer>(m_memoryManager, m_width, m_height, false, m_samples);
    m_depthFramebuffer->setDepth(&m_depthImageView[0], VK_FORMAT_D32_SFLOAT);
    m_depthFramebuffer->initializeFramebuffer(&m_depthRenderPass);
}

void
VulkanRenderer::deleteFramebuffers()
{
    // The framebuffers/command buffers may still be in use
    vkDeviceWaitIdle(m_renderDevice);

    // Depth buffer
    for (uint32_t i = 0; i < m_mipLevels; i++)
    {
        vkDestroyImageView(m_renderDevice, m_depthImageView[i], nullptr);
    }

    // HDR buffers
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < m_HDRImageView[i].size(); j++)
        {
            vkDestroyImageView(m_renderDevice, m_HDRImageView[i][j], nullptr);
        }
    }

    // Normal buffer
    vkDestroyImageView(m_renderDevice, m_normalImageView, nullptr);

    // AO buffer
    vkDestroyImageView(m_renderDevice, m_halfAOImageView[0], nullptr);
    vkDestroyImageView(m_renderDevice, m_halfAOImageView[1], nullptr);

    for (auto imageView : m_swapchainImageViews)
    {
        vkDestroyImageView(m_renderDevice, imageView, nullptr);
    }

    // Delete all post processing resources
    for (auto postProcess : m_postProcessingChain->m_postProcesses)
    {
        postProcess->m_framebuffer->clear(&m_renderDevice);
    }

    // Delete all HDR resources
    for (auto pass : m_HDRTonemaps)
    {
        pass->m_framebuffer->clear(&m_renderDevice);
    }

    // Delete all AO resources
    for (auto pass : m_ssao)
    {
        pass->m_framebuffer->clear(&m_renderDevice);
    }

    // Delete all drawing resources
    vkDestroyFramebuffer(m_renderDevice, m_opaqueFramebuffer->m_framebuffer, nullptr);
    vkDestroyFramebuffer(m_renderDevice, m_depthFramebuffer->m_framebuffer, nullptr);
    vkDestroyFramebuffer(m_renderDevice, m_decalFramebuffer->m_framebuffer, nullptr);
}

void
VulkanRenderer::renderFrame()
{
    m_frameNumber++;

    // The swapchain contains multiple buffers, so get one that is available (i.e., not currently being written to)
    uint32_t nextImageIndex;
    vkAcquireNextImageKHR(m_renderDevice, *m_swapchain, UINT64_MAX, m_readyToRender, VK_NULL_HANDLE, &nextImageIndex);

    this->loadAllGeometry();

    // Update global uniforms
    this->updateGlobalUniforms(nextImageIndex);

    // Update local uniforms
    for (unsigned int renderDelegateIndex = 0; renderDelegateIndex < m_renderDelegates.size(); renderDelegateIndex++)
    {
        if (m_renderDelegates[renderDelegateIndex]->getGeometry()->getType() == Geometry::Type::DecalPool)
        {
            auto decalPool = std::dynamic_pointer_cast<VulkanDecalRenderDelegate>(m_renderDelegates[renderDelegateIndex]);
            decalPool->update(nextImageIndex, m_scene->getCamera());
        }
        m_renderDelegates[renderDelegateIndex]->update(nextImageIndex);
    }

    // Wait until command buffer is done so that we can write to it again
    vkWaitForFences(m_renderDevice, 1, &m_commandBufferSubmit[nextImageIndex], VK_TRUE, UINT64_MAX);
    vkResetFences(m_renderDevice, 1, &m_commandBufferSubmit[nextImageIndex]);

    VkCommandBufferBeginInfo commandBufferBeginInfo;
    commandBufferBeginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    commandBufferBeginInfo.pNext = nullptr;
    commandBufferBeginInfo.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    commandBufferBeginInfo.pInheritanceInfo = nullptr;

    vkBeginCommandBuffer(m_renderCommandBuffer[nextImageIndex], &commandBufferBeginInfo);

    VkRect2D renderArea;
    renderArea.offset = { 0, 0 };
    renderArea.extent = { m_width, m_height };

    std::array<VkClearValue, 4> clearValues;
    clearValues[0].color = { { (float)m_backgroundColor[0], (float)m_backgroundColor[1], (float)m_backgroundColor[2], 1 } }; // Color
    clearValues[1].depthStencil = { { 1.0 }, { 0 } }; // Depth
    clearValues[2].color = { { 0, 0, 0, 0 } }; // Normal
    clearValues[3].color = { { 0, 0, 0, 0 } }; // Specular

    // Do buffer transfers
    for (unsigned int renderDelegateIndex = 0; renderDelegateIndex < m_renderDelegates.size(); renderDelegateIndex++)
    {
        auto buffers = m_renderDelegates[renderDelegateIndex]->getBuffer().get();
        buffers->uploadBuffers(m_renderCommandBuffer[nextImageIndex]);
    }

    VkDeviceSize deviceSize = { 0 };

    // Pass 0: Render opaque shadows
    for (size_t i = 0; i < m_shadowPasses.size(); i++)
    {
        VkRect2D shadowRenderArea = { {0, 0}, {m_shadowMapResolution, m_shadowMapResolution} };

        VkRenderPassBeginInfo shadowPassBeginInfo;
        shadowPassBeginInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
        shadowPassBeginInfo.pNext = nullptr;
        shadowPassBeginInfo.renderPass = m_shadowPasses[i];
        shadowPassBeginInfo.framebuffer = m_shadowFramebuffers[i]->m_framebuffer;
        shadowPassBeginInfo.renderArea = shadowRenderArea;
        shadowPassBeginInfo.clearValueCount = 1;
        shadowPassBeginInfo.pClearValues = &clearValues[1]; // depth buffer

        vkCmdBeginRenderPass(m_renderCommandBuffer[nextImageIndex], &shadowPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

        for (unsigned int renderDelegateIndex = 0; renderDelegateIndex < m_renderDelegates.size(); renderDelegateIndex++)
        {
            if (m_renderDelegates[renderDelegateIndex]->getGeometry()->getType() == Geometry::Type::DecalPool)
            {
                continue;
            }

            auto material = m_renderDelegates[renderDelegateIndex]->m_shadowMaterial;

            if (!m_renderDelegates[renderDelegateIndex]->getGeometry()->getRenderMaterial()->getCastsShadows())
            {
                continue;
            }

            vkCmdBindPipeline(m_renderCommandBuffer[nextImageIndex], VK_PIPELINE_BIND_POINT_GRAPHICS, material->m_pipeline);

            vkCmdBindDescriptorSets(m_renderCommandBuffer[nextImageIndex],
                VK_PIPELINE_BIND_POINT_GRAPHICS,
                material->m_pipelineLayout, 0, (uint32_t)material->m_descriptorSets.size(),
                &material->m_descriptorSets[0], 0, &m_dynamicOffsets);

            auto buffers = m_renderDelegates[renderDelegateIndex]->getBuffer().get();

            vkCmdPushConstants(m_renderCommandBuffer[nextImageIndex], material->m_pipelineLayout, VK_SHADER_STAGE_VERTEX_BIT, 0, 64, (void*)&m_lightMatrices[i]);
            buffers->bindBuffers(&m_renderCommandBuffer[nextImageIndex], nextImageIndex);
            vkCmdDrawIndexed(m_renderCommandBuffer[nextImageIndex], buffers->m_numIndices, 1, 0, 0, 0);
        }
        vkCmdEndRenderPass(m_renderCommandBuffer[nextImageIndex]);
        VulkanAttachmentBarriers::addDepthAttachmentBarrier(&m_renderCommandBuffer[nextImageIndex],
                        m_renderQueueFamily,
                        m_shadowMaps->getImage());
    }

    // Pass 1: Depth pre-pass
    VkRenderPassBeginInfo depthRenderPassBeginInfo;
    depthRenderPassBeginInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
    depthRenderPassBeginInfo.pNext = nullptr;
    depthRenderPassBeginInfo.renderPass = m_depthRenderPass;
    depthRenderPassBeginInfo.framebuffer = m_depthFramebuffer->m_framebuffer;
    depthRenderPassBeginInfo.renderArea = renderArea;
    depthRenderPassBeginInfo.clearValueCount = 1;
    depthRenderPassBeginInfo.pClearValues = &clearValues[1];

    vkCmdBeginRenderPass(m_renderCommandBuffer[nextImageIndex], &depthRenderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

    for (unsigned int renderDelegateIndex = 0; renderDelegateIndex < m_renderDelegates.size(); renderDelegateIndex++)
    {
        if (m_renderDelegates[renderDelegateIndex]->getGeometry()->getType() == Geometry::Type::DecalPool)
        {
            continue;
        }

        auto material = m_renderDelegates[renderDelegateIndex]->m_depthMaterial;
        vkCmdBindPipeline(m_renderCommandBuffer[nextImageIndex], VK_PIPELINE_BIND_POINT_GRAPHICS, material->m_pipeline);
        this->setCommandBufferState(&m_renderCommandBuffer[nextImageIndex], m_width, m_height);

        vkCmdBindDescriptorSets(m_renderCommandBuffer[nextImageIndex],
            VK_PIPELINE_BIND_POINT_GRAPHICS,
            material->m_pipelineLayout, 0, (uint32_t)material->m_descriptorSets.size(),
            &material->m_descriptorSets[0], 0, &m_dynamicOffsets);

        auto buffers = m_renderDelegates[renderDelegateIndex]->getBuffer().get();

        buffers->bindBuffers(&m_renderCommandBuffer[nextImageIndex], nextImageIndex);
        vkCmdDrawIndexed(m_renderCommandBuffer[nextImageIndex], buffers->m_numIndices, 1, 0, 0, 0);
    }
    vkCmdEndRenderPass(m_renderCommandBuffer[nextImageIndex]);
    VulkanAttachmentBarriers::addDepthAttachmentBarrier(&m_renderCommandBuffer[nextImageIndex], m_renderQueueFamily, m_depthImage[0]->getImage());

    // Pass 1 - 4: AO processing
    VkRenderPassBeginInfo aoRenderPassBeginInfo;
    aoRenderPassBeginInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
    aoRenderPassBeginInfo.pNext = nullptr;
    aoRenderPassBeginInfo.clearValueCount = (uint32_t)clearValues.size();
    aoRenderPassBeginInfo.pClearValues = &clearValues[0];

    for (unsigned int postProcessIndex = 0; postProcessIndex < m_ssao.size(); postProcessIndex++)
    {
        auto postProcess = m_ssao[postProcessIndex];

        auto framebuffer = postProcess->m_framebuffer;
        aoRenderPassBeginInfo.renderPass =
            postProcess->m_renderPass;
        aoRenderPassBeginInfo.framebuffer = framebuffer->m_framebuffer;
        aoRenderPassBeginInfo.clearValueCount = (uint32_t)framebuffer->m_attachments.size();
        aoRenderPassBeginInfo.renderArea.offset = { 0, 0 };
        aoRenderPassBeginInfo.renderArea.extent = { framebuffer->m_width, framebuffer->m_height };

        vkCmdBeginRenderPass(m_renderCommandBuffer[nextImageIndex], &aoRenderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

        vkCmdPushConstants(m_renderCommandBuffer[nextImageIndex], postProcess->m_pipelineLayout,
            VK_SHADER_STAGE_FRAGMENT_BIT, 0, 128, (void*)postProcess->m_pushConstantData);

        vkCmdBindPipeline(m_renderCommandBuffer[nextImageIndex], VK_PIPELINE_BIND_POINT_GRAPHICS,
            postProcess->m_pipeline);
        this->setCommandBufferState(&m_renderCommandBuffer[nextImageIndex], framebuffer->m_width, framebuffer->m_height);

        vkCmdBindDescriptorSets(m_renderCommandBuffer[nextImageIndex],
            VK_PIPELINE_BIND_POINT_GRAPHICS,
            postProcess->m_pipelineLayout, 0,
            (uint32_t)postProcess->m_descriptorSets.size(),
            &postProcess->m_descriptorSets[0], 0, &m_dynamicOffsets);

        auto buffers = postProcess->m_vertexBuffer;
        buffers->bindBuffers(&m_renderCommandBuffer[nextImageIndex], 0);
        vkCmdDrawIndexed(m_renderCommandBuffer[nextImageIndex], buffers->m_numIndices, 1, 0, 0, 0);

        vkCmdEndRenderPass(m_renderCommandBuffer[nextImageIndex]);
    }

    // Pass 2: Render opaque geometry
    VkRenderPassBeginInfo opaqueRenderPassBeginInfo;
    opaqueRenderPassBeginInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
    opaqueRenderPassBeginInfo.pNext = nullptr;
    opaqueRenderPassBeginInfo.renderPass = m_opaqueRenderPass;
    opaqueRenderPassBeginInfo.framebuffer = m_opaqueFramebuffer->m_framebuffer;
    opaqueRenderPassBeginInfo.renderArea = renderArea;
    opaqueRenderPassBeginInfo.clearValueCount = (uint32_t)clearValues.size();
    opaqueRenderPassBeginInfo.pClearValues = &clearValues[0];

    vkCmdBeginRenderPass(m_renderCommandBuffer[nextImageIndex], &opaqueRenderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

    for (unsigned int renderDelegateIndex = 0; renderDelegateIndex < m_renderDelegates.size(); renderDelegateIndex++)
    {
        if (m_renderDelegates[renderDelegateIndex]->getGeometry()->getType() == Geometry::Type::DecalPool)
        {
            continue;
        }

        auto material = m_renderDelegates[renderDelegateIndex]->m_material;
        vkCmdBindPipeline(m_renderCommandBuffer[nextImageIndex], VK_PIPELINE_BIND_POINT_GRAPHICS, material->m_pipeline);
        this->setCommandBufferState(&m_renderCommandBuffer[nextImageIndex], m_width, m_height);

        vkCmdBindDescriptorSets(m_renderCommandBuffer[nextImageIndex],
            VK_PIPELINE_BIND_POINT_GRAPHICS,
            material->m_pipelineLayout, 0, (uint32_t)material->m_descriptorSets.size(),
            &material->m_descriptorSets[0], 0, &m_dynamicOffsets);

        auto buffers = m_renderDelegates[renderDelegateIndex]->getBuffer().get();

        buffers->bindBuffers(&m_renderCommandBuffer[nextImageIndex], nextImageIndex);
        vkCmdDrawIndexed(m_renderCommandBuffer[nextImageIndex], buffers->m_numIndices, 1, 0, 0, 0);
    }
    vkCmdEndRenderPass(m_renderCommandBuffer[nextImageIndex]);

    // Pass 2: Render decals
    VkRenderPassBeginInfo decalRenderPassBeginInfo;
    decalRenderPassBeginInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
    decalRenderPassBeginInfo.pNext = nullptr;
    decalRenderPassBeginInfo.renderPass = m_decalRenderPass;
    decalRenderPassBeginInfo.framebuffer = m_decalFramebuffer->m_framebuffer;
    decalRenderPassBeginInfo.renderArea = renderArea;
    decalRenderPassBeginInfo.clearValueCount = 0;
    decalRenderPassBeginInfo.pClearValues = &clearValues[0];
    vkCmdBeginRenderPass(m_renderCommandBuffer[nextImageIndex], &decalRenderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

    for (unsigned int renderDelegateIndex = 0; renderDelegateIndex < m_renderDelegates.size(); renderDelegateIndex++)
    {
        if (m_renderDelegates[renderDelegateIndex]->getGeometry()->getType() != Geometry::Type::DecalPool)
        {
            continue;
        }

        auto geometry = std::dynamic_pointer_cast<DecalPool>(m_renderDelegates[renderDelegateIndex]->getGeometry());
        auto material = m_renderDelegates[renderDelegateIndex]->m_material;
        vkCmdBindPipeline(m_renderCommandBuffer[nextImageIndex], VK_PIPELINE_BIND_POINT_GRAPHICS, material->m_pipeline);
        this->setCommandBufferState(&m_renderCommandBuffer[nextImageIndex], m_width, m_height);

        vkCmdBindDescriptorSets(m_renderCommandBuffer[nextImageIndex],
            VK_PIPELINE_BIND_POINT_GRAPHICS,
            material->m_pipelineLayout, 0, (uint32_t)material->m_descriptorSets.size(),
            &material->m_descriptorSets[0], 0, &m_dynamicOffsets);

        auto buffers = m_renderDelegates[renderDelegateIndex]->getBuffer().get();

        buffers->bindBuffers(&m_renderCommandBuffer[nextImageIndex], nextImageIndex);
        vkCmdDrawIndexed(m_renderCommandBuffer[nextImageIndex], buffers->m_numIndices, geometry->getNumDecals(), 0, 0, 0);
    }

    vkCmdEndRenderPass(m_renderCommandBuffer[nextImageIndex]);
    vkEndCommandBuffer(m_renderCommandBuffer[nextImageIndex]);

    vkBeginCommandBuffer(m_postProcessingCommandBuffer[nextImageIndex], &commandBufferBeginInfo);

    // Pass 3 to N - 1: Post processing
    for (unsigned int postProcessIndex = 0; postProcessIndex < m_postProcessingChain->m_postProcesses.size(); postProcessIndex++)
    {
        clearValues[0].color = { { 1.0, 0.0, 0.0, 1 } }; // Color

        auto postProcess = m_postProcessingChain->m_postProcesses[postProcessIndex];

        auto postProcessRenderPassBeginInfo = opaqueRenderPassBeginInfo;
        auto framebuffer = postProcess->m_framebuffer;
        postProcessRenderPassBeginInfo.renderPass =
            postProcess->m_renderPass;
        postProcessRenderPassBeginInfo.framebuffer = framebuffer->m_framebuffer;
        postProcessRenderPassBeginInfo.clearValueCount = (uint32_t)framebuffer->m_attachments.size();
        postProcessRenderPassBeginInfo.renderArea.extent = { framebuffer->m_width, framebuffer->m_height };

        vkCmdBeginRenderPass(m_postProcessingCommandBuffer[nextImageIndex], &postProcessRenderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

        vkCmdPushConstants(m_postProcessingCommandBuffer[nextImageIndex], postProcess->m_pipelineLayout,
            VK_SHADER_STAGE_FRAGMENT_BIT, 0, 128, (void*)postProcess->m_pushConstantData);

        vkCmdBindPipeline(m_postProcessingCommandBuffer[nextImageIndex], VK_PIPELINE_BIND_POINT_GRAPHICS,
            postProcess->m_pipeline);
        this->setCommandBufferState(&m_postProcessingCommandBuffer[nextImageIndex], framebuffer->m_width, framebuffer->m_height);

        vkCmdBindDescriptorSets(m_postProcessingCommandBuffer[nextImageIndex],
            VK_PIPELINE_BIND_POINT_GRAPHICS,
            postProcess->m_pipelineLayout, 0,
            (uint32_t)postProcess->m_descriptorSets.size(),
            &postProcess->m_descriptorSets[0], 0, &m_dynamicOffsets);

        auto buffers = postProcess->m_vertexBuffer;
        buffers->bindBuffers(&m_postProcessingCommandBuffer[nextImageIndex], 0);
        vkCmdDrawIndexed(m_postProcessingCommandBuffer[nextImageIndex], buffers->m_numIndices, 1, 0, 0, 0);

        vkCmdEndRenderPass(m_postProcessingCommandBuffer[nextImageIndex]);
    }

    // Pass N: HDR tonemap (this is special because of the swapchain)
    {
        auto postProcessRenderPassBeginInfo = opaqueRenderPassBeginInfo;
        postProcessRenderPassBeginInfo.renderPass = m_HDRTonemaps[nextImageIndex]->m_renderPass;
        postProcessRenderPassBeginInfo.framebuffer = m_HDRTonemaps[nextImageIndex]->m_framebuffer->m_framebuffer;
        postProcessRenderPassBeginInfo.clearValueCount = (uint32_t)m_HDRTonemaps[nextImageIndex]->m_framebuffer->m_attachments.size();

        vkCmdBeginRenderPass(m_postProcessingCommandBuffer[nextImageIndex], &postProcessRenderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

        vkCmdBindPipeline(m_postProcessingCommandBuffer[nextImageIndex], VK_PIPELINE_BIND_POINT_GRAPHICS, m_HDRTonemaps[nextImageIndex]->m_pipeline);
        this->setCommandBufferState(&m_postProcessingCommandBuffer[nextImageIndex],
            m_HDRTonemaps[nextImageIndex]->m_framebuffer->m_width,
            m_HDRTonemaps[nextImageIndex]->m_framebuffer->m_height);

        vkCmdBindDescriptorSets(m_postProcessingCommandBuffer[nextImageIndex],
            VK_PIPELINE_BIND_POINT_GRAPHICS,
            m_HDRTonemaps[nextImageIndex]->m_pipelineLayout, 0, (uint32_t)m_HDRTonemaps[nextImageIndex]->m_descriptorSets.size(),
            &m_HDRTonemaps[nextImageIndex]->m_descriptorSets[0], 0, &m_dynamicOffsets);

        auto buffers = m_HDRTonemaps[nextImageIndex]->m_vertexBuffer;
        buffers->bindBuffers(&m_postProcessingCommandBuffer[nextImageIndex], 0);
        vkCmdDrawIndexed(m_postProcessingCommandBuffer[nextImageIndex], buffers->m_numIndices, 1, 0, 0, 0);

        vkCmdEndRenderPass(m_postProcessingCommandBuffer[nextImageIndex]);
    }

    vkEndCommandBuffer(m_postProcessingCommandBuffer[nextImageIndex]);

    VkCommandBuffer commandBuffers[2];
    commandBuffers[0] = m_renderCommandBuffer[nextImageIndex];
    commandBuffers[1] = m_postProcessingCommandBuffer[nextImageIndex];

    VkPipelineStageFlags stageWaitFlags = VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT;
    VkSubmitInfo submitInfo[2];
    submitInfo[0].sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
    submitInfo[0].pNext = nullptr;
    submitInfo[0].waitSemaphoreCount = 1;
    submitInfo[0].pWaitSemaphores = &m_readyToRender;
    submitInfo[0].pWaitDstStageMask = &stageWaitFlags;
    submitInfo[0].commandBufferCount = 1;
    submitInfo[0].pCommandBuffers = &commandBuffers[0];
    submitInfo[0].signalSemaphoreCount = 1;
    submitInfo[0].pSignalSemaphores = &m_drawingComplete;

    submitInfo[1] = submitInfo[0];
    submitInfo[1].pWaitSemaphores = &m_drawingComplete;
    submitInfo[1].pWaitDstStageMask = &stageWaitFlags;
    submitInfo[1].pCommandBuffers = &commandBuffers[1];
    submitInfo[1].pSignalSemaphores = &m_presentImages;

    // Submit command buffers
    vkQueueSubmit(m_renderQueue, 2, submitInfo, m_commandBufferSubmit[nextImageIndex]);

    VkSwapchainKHR swapchains[] = {*m_swapchain};

    VkPresentInfoKHR presentInfo;
    presentInfo.sType = VK_STRUCTURE_TYPE_PRESENT_INFO_KHR;
    presentInfo.pNext = nullptr;
    presentInfo.waitSemaphoreCount = 1;
    presentInfo.pWaitSemaphores = &m_presentImages;
    presentInfo.swapchainCount = 1;
    presentInfo.pSwapchains = swapchains;
    presentInfo.pImageIndices = &nextImageIndex; // change this in the future for multi-screen rendering
    presentInfo.pResults = nullptr;

    // Display backbuffer
    vkQueuePresentKHR(m_renderQueue, &presentInfo);
}

void
VulkanRenderer::setupSynchronization()
{
    VkSemaphoreCreateInfo semaphoreInfo;
    semaphoreInfo.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;
    semaphoreInfo.pNext = nullptr;
    semaphoreInfo.flags = 0;

    vkCreateSemaphore(m_renderDevice, &semaphoreInfo, nullptr, &m_readyToRender);
    vkCreateSemaphore(m_renderDevice, &semaphoreInfo, nullptr, &m_presentImages);
    vkCreateSemaphore(m_renderDevice, &semaphoreInfo, nullptr, &m_drawingComplete);

    VkFenceCreateInfo fenceInfo;
    fenceInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
    fenceInfo.pNext = nullptr;
    fenceInfo.flags = VK_FENCE_CREATE_SIGNALED_BIT;

    m_commandBufferSubmit.resize(m_buffering);
    for (uint32_t i = 0; i < m_buffering; i++)
    {
        vkCreateFence(m_renderDevice, &fenceInfo, nullptr, &m_commandBufferSubmit[i]);
    }
}

void
VulkanRenderer::loadAllGeometry()
{
    // Add new objects
    for (auto sceneObject : m_scene->getSceneObjects())
    {
        auto type = sceneObject->getType();
        auto geometry = sceneObject->getVisualGeometry();

        if (geometry && !geometry->m_renderDelegateCreated)
        {
            auto renderDelegate = this->loadGeometry(geometry, type);
        }
    }
}

std::shared_ptr<VulkanRenderDelegate>
VulkanRenderer::loadGeometry(std::shared_ptr<Geometry> geometry, SceneObject::Type type)
{
    auto renderDelegate = VulkanRenderDelegate::make_delegate(geometry, type, m_memoryManager);
    if (renderDelegate != nullptr)
    {
        m_renderDelegates.push_back(renderDelegate);
        renderDelegate->getBuffer()->initializeBuffers(m_memoryManager);
        renderDelegate->m_material->initialize(this);

        if (!renderDelegate->getGeometry()->getRenderMaterial()->isDecal())
        {
            renderDelegate->m_shadowMaterial->initialize(this);
            renderDelegate->m_depthMaterial->initialize(this);
        }
    }
    return renderDelegate;
}

void
VulkanRenderer::setupMemoryManager()
{
    m_memoryManager.setup(&m_renderPhysicalDevice);
    m_memoryManager.m_device = m_renderDevice;
    m_memoryManager.m_queueFamilyIndex = m_renderQueueFamily;
    m_memoryManager.m_transferCommandBuffer = &m_renderCommandBuffer[0];
    m_memoryManager.m_transferQueue = &m_renderQueue;
}

void
VulkanRenderer::createGlobalUniformBuffers()
{
    m_globalVertexUniformBuffer = std::make_shared<VulkanUniformBuffer>(m_memoryManager, (uint32_t)sizeof(VulkanGlobalVertexUniforms));
    m_globalFragmentUniformBuffer = std::make_shared<VulkanUniformBuffer>(m_memoryManager, (uint32_t)sizeof(VulkanGlobalFragmentUniforms));
}

void
VulkanRenderer::initializePostProcesses()
{
    std::vector<VkPipeline> graphicsPipelines;
    std::vector<VkGraphicsPipelineCreateInfo> graphicsPipelinesInfo;
    m_postProcessingChain = std::make_shared<VulkanPostProcessingChain>(this);

    // HDR pipeline creation
    for (uint32_t i = 0; i < m_swapchainImageCount; i++)
    {
        m_HDRTonemaps[i] = std::make_shared<VulkanPostProcess>(this, 0, true);
        m_HDRTonemaps[i]->addInputImage(&m_swapchainImageSamplers[i], &m_HDRImageView[m_postProcessingChain->m_lastOutput][0]);
        m_HDRTonemaps[i]->m_framebuffer->setColor(&m_swapchainImageViews[i], VK_FORMAT_B8G8R8A8_SRGB);
        m_HDRTonemaps[i]->initialize(this, "./Shaders/VulkanShaders/PostProcessing/HDR_tonemap_frag.spv");

        graphicsPipelines.push_back(m_HDRTonemaps[i]->m_pipeline);
        graphicsPipelinesInfo.push_back(m_HDRTonemaps[i]->m_graphicsPipelineInfo);
    }

    // AO pipeline creation
    m_ssao.resize(4);

    // Textures
    if (!m_noiseTexture)
    {
        m_noiseTexture = std::make_shared<Texture>("noise", Texture::Type::NONE);
        m_noiseTextureDelegate = std::make_shared<VulkanTextureDelegate>(m_memoryManager, m_noiseTexture);
    }

    m_ssao[0] = std::make_shared<VulkanPostProcess>(this, 1);
    m_ssao[0]->addInputImage(&m_HDRImageSampler, &m_depthImageView[0], VK_IMAGE_LAYOUT_DEPTH_STENCIL_READ_ONLY_OPTIMAL);
    m_ssao[0]->m_framebuffer->setColor(&m_depthImageView[1], VK_FORMAT_R32_SFLOAT);
    m_ssao[0]->initialize(this, "./Shaders/VulkanShaders/PostProcessing/depth_downscale_frag.spv");

    m_ssao[1] = std::make_shared<VulkanPostProcess>(this, 1);
    m_ssao[1]->addInputImage(&m_HDRImageSampler, &m_depthImageView[1]);
    m_ssao[1]->addInputImage(&m_noiseTextureDelegate->m_sampler, &m_noiseTextureDelegate->m_imageView);
    m_ssao[1]->m_framebuffer->setColor(&m_halfAOImageView[0], VK_FORMAT_R8_UNORM);
    m_ssao[1]->m_pushConstantData[0] = m_fov;
    m_ssao[1]->m_pushConstantData[1] = 0.1;
    m_ssao[1]->m_pushConstantData[2] = m_nearPlane;
    m_ssao[1]->m_pushConstantData[3] = m_farPlane;
    m_ssao[1]->m_pushConstantData[4] = 6;
    m_ssao[1]->m_pushConstantData[5] = m_width / 2;
    m_ssao[1]->m_pushConstantData[6] = m_height / 2;
    m_ssao[1]->initialize(this, "./Shaders/VulkanShaders/PostProcessing/ao_frag.spv");

    m_ssao[2] = std::make_shared<VulkanPostProcess>(this, 1);
    m_ssao[2]->addInputImage(&m_HDRImageSampler, &m_halfAOImageView[0]);
    m_ssao[2]->addInputImage(&m_HDRImageSampler, &m_depthImageView[1]);
    m_ssao[2]->m_framebuffer->setColor(&m_halfAOImageView[1], VK_FORMAT_R8_UNORM);
    m_ssao[2]->m_pushConstantData[0] = std::max(m_width >> 1, 1u);
    m_ssao[2]->m_pushConstantData[1] = std::max(m_height >> 1, 1u);
    m_ssao[2]->m_pushConstantData[2] = m_nearPlane;
    m_ssao[2]->m_pushConstantData[3] = m_farPlane;
    m_ssao[2]->m_pushConstantData[4] = 2;
    VulkanPostProcessingChain::calculateBlurValuesLinear(2,
        &m_ssao[2]->m_pushConstantData[5],
        &m_ssao[2]->m_pushConstantData[10]);
    m_ssao[2]->initialize(this, "./Shaders/VulkanShaders/PostProcessing/bilateral_blur_horizontal_frag.spv");

    m_ssao[3] = std::make_shared<VulkanPostProcess>(this, 1);
    m_ssao[3]->addInputImage(&m_HDRImageSampler, &m_halfAOImageView[1]);
    m_ssao[3]->addInputImage(&m_HDRImageSampler, &m_depthImageView[1]);
    m_ssao[3]->m_framebuffer->setColor(&m_halfAOImageView[0], VK_FORMAT_R8_UNORM);
    m_ssao[3]->m_pushConstantData[0] = std::max(m_width >> 1, 1u);
    m_ssao[3]->m_pushConstantData[1] = std::max(m_height >> 1, 1u);
    m_ssao[3]->m_pushConstantData[2] = m_nearPlane;
    m_ssao[3]->m_pushConstantData[3] = m_farPlane;
    m_ssao[3]->m_pushConstantData[4] = 2;
    VulkanPostProcessingChain::calculateBlurValuesLinear(2,
        &m_ssao[3]->m_pushConstantData[5],
        &m_ssao[3]->m_pushConstantData[10]);
    m_ssao[3]->initialize(this, "./Shaders/VulkanShaders/PostProcessing/bilateral_blur_vertical_frag.spv");

    for (int i = 0; i < 4; i++)
    {
        graphicsPipelines.push_back(m_ssao[i]->m_pipeline);
        graphicsPipelinesInfo.push_back(m_ssao[i]->m_graphicsPipelineInfo);
    }

    // Post processing pipeline creation
    for (int i = 0; i < m_postProcessingChain->m_postProcesses.size(); i++)
    {
        graphicsPipelines.push_back(m_postProcessingChain->m_postProcesses[i]->m_pipeline);
        graphicsPipelinesInfo.push_back(m_postProcessingChain->m_postProcesses[i]->m_graphicsPipelineInfo);
    }

    vkCreateGraphicsPipelines(m_renderDevice,
        m_pipelineCache,
        (uint32_t)graphicsPipelines.size(),
        &graphicsPipelinesInfo[0],
        nullptr,
        &graphicsPipelines[0]);

    int index = 0;
    for (uint32_t i = 0; i < m_swapchainImageCount; i++)
    {
        m_HDRTonemaps[i]->m_pipeline = graphicsPipelines[index];
        index++;
    }

    for (uint32_t i = 0; i < 4; i++)
    {
        m_ssao[i]->m_pipeline = graphicsPipelines[index];
        index++;
    }

    for (int i = 0; i < m_postProcessingChain->m_postProcesses.size(); i++)
    {
        m_postProcessingChain->m_postProcesses[i]->m_pipeline = graphicsPipelines[index];
        index++;
    }
}

void
VulkanRenderer::updateGlobalUniforms(uint32_t frameIndex)
{
    // Vertex uniforms
    {
        // Projection matrix
        auto camera = m_scene->getCamera();
        m_fov = (float)glm::radians(camera->getViewAngle());
        m_globalVertexUniforms.projectionMatrix = glm::perspective(m_fov, (float)(m_width) / (float)(m_height), m_nearPlane, m_farPlane);
        glm::mat4 correctionMatrix; // for Vulkan rendering
        correctionMatrix[1][1] = -1;
        correctionMatrix[2][2] = 0.5;
        correctionMatrix[3][2] = 0.5;
        m_globalVertexUniforms.projectionMatrix *= correctionMatrix;

        // View matrix
        auto eye = glm::tvec3<float>(camera->getPosition().x(), camera->getPosition().y(), camera->getPosition().z());
        auto center = glm::tvec3<float>(camera->getFocalPoint().x(), camera->getFocalPoint().y(), camera->getFocalPoint().z());
        auto up = glm::tvec3<float>(camera->getViewUp().x(), camera->getViewUp().y(), camera->getViewUp().z());
        m_globalVertexUniforms.cameraPosition = glm::vec4(camera->getPosition().x(), camera->getPosition().y(), camera->getPosition().z(), 0.0);
        m_globalVertexUniforms.viewMatrix = glm::lookAt(eye, center, up);
    }

    // Lights uniforms
    {
        // Update lights - need to optimize
        auto lights = m_scene->getLights();
        for (int i = 0; i < lights.size(); i++)
        {
            // Only supports directional lights right now
            auto focalPoint = lights[i]->getFocalPoint();
            auto position = Vec3d(0,0,0);
            int type = 1;
            int shadowMapIndex = -1;

            if (lights[i]->getType() == LightType::POINT_LIGHT || lights[i]->getType() == LightType::SPOT_LIGHT)
            {
                position = std::static_pointer_cast<PointLight>(lights[i])->getPosition();
                type = 2;
            }

            m_globalFragmentUniforms.lights[i].position = glm::vec4(position.x(), position.y(), position.z(), 1);

            m_globalFragmentUniforms.lights[i].direction.x = focalPoint.x() - position.x();
            m_globalFragmentUniforms.lights[i].direction.y = focalPoint.y() - position.y();
            m_globalFragmentUniforms.lights[i].direction.z = focalPoint.z() - position.z();

            m_globalFragmentUniforms.lights[i].direction = glm::normalize(m_globalFragmentUniforms.lights[i].direction);

            Color lightColor = lights[i]->getColor();
            m_globalFragmentUniforms.lights[i].color = glm::vec4(lightColor.r, lightColor.g, lightColor.b, 1.0);

            if (lights[i]->getType() == LightType::SPOT_LIGHT)
            {
                m_globalFragmentUniforms.lights[i].direction.a =
                    glm::radians(std::static_pointer_cast<SpotLight>(lights[i])->getSpotAngle());
                type = 3;
            }

            if (lights[i]->getType() == LightType::DIRECTIONAL_LIGHT)
            {
                shadowMapIndex = std::static_pointer_cast<DirectionalLight>(lights[i])->m_shadowMapIndex;
            }

            m_globalFragmentUniforms.lights[i].color.a = lights[i]->getIntensity();

            m_globalFragmentUniforms.lights[i].state.x = type;
            m_globalFragmentUniforms.lights[i].state.y = shadowMapIndex;
        }

        memcpy(&m_globalVertexUniforms.lights, &m_globalFragmentUniforms.lights, sizeof(m_globalFragmentUniforms.lights));

        m_globalFragmentUniforms.inverseViewMatrix = glm::inverse(m_globalVertexUniforms.viewMatrix);
        m_globalFragmentUniforms.inverseProjectionMatrix = glm::inverse(m_globalVertexUniforms.projectionMatrix);
        m_globalFragmentUniforms.resolution = glm::vec4(m_width, m_height, m_shadowMapResolution, 0);

        for (size_t i = 0; i < m_shadowLights.size(); i++)
        {
            auto light = m_shadowLights[i];
            auto camera = m_scene->getCamera();

            auto shadowRange = light->m_shadowRange;
            auto shadowCenter = light->m_shadowCenter;

            m_lightMatrices[i] = glm::ortho(-shadowRange, shadowRange, -shadowRange, shadowRange, -shadowRange, shadowRange);
            glm::mat4 correctionMatrix; // for Vulkan rendering
            correctionMatrix[1][1] = -1;
            correctionMatrix[2][2] = 0.5;
            correctionMatrix[3][2] = 0.5;
            m_lightMatrices[i] *= correctionMatrix;
            auto eye = glm::tvec3<float>(shadowCenter.x(), shadowCenter.y(), shadowCenter.z());
            auto center = glm::tvec3<float>(light->getFocalPoint().x(), light->getFocalPoint().y(), light->getFocalPoint().z()) + eye;
            auto offset = glm::normalize(eye - center) * shadowRange;
            center += offset;
            eye += offset;
            auto up = glm::tvec3<float>(0, 1, 0);
            m_lightMatrices[i] *= glm::lookAt(eye, center, up);
            m_globalFragmentUniforms.lightMatrices[i] = m_lightMatrices[i];
        }
    }

    m_globalVertexUniformBuffer->updateUniforms(sizeof(VulkanGlobalVertexUniforms), &m_globalVertexUniforms, frameIndex);
    m_globalFragmentUniformBuffer->updateUniforms(sizeof(VulkanGlobalFragmentUniforms), &m_globalFragmentUniforms, frameIndex);
}

void
VulkanRenderer::createShadowMaps(uint32_t resolution)
{
    // shadow maps image
    uint32_t numShadows = 0;
    for (auto light : m_scene->getLights())
    {
        if (light->getType() == LightType::DIRECTIONAL_LIGHT)
        {
            numShadows++;
        }
    }

    VkImageCreateInfo shadowMapsInfo;
    shadowMapsInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
    shadowMapsInfo.pNext = nullptr;
    shadowMapsInfo.flags = 0;
    shadowMapsInfo.imageType = VK_IMAGE_TYPE_2D;
    shadowMapsInfo.format = VK_FORMAT_D32_SFLOAT;
    shadowMapsInfo.extent = { resolution, resolution, 1 };
    shadowMapsInfo.mipLevels = 1;
    shadowMapsInfo.arrayLayers = std::max(numShadows, (uint32_t)1);
    shadowMapsInfo.samples = m_samples;
    shadowMapsInfo.tiling = VK_IMAGE_TILING_OPTIMAL;
    shadowMapsInfo.usage = VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT
                           | VK_IMAGE_USAGE_SAMPLED_BIT;
    shadowMapsInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    shadowMapsInfo.queueFamilyIndexCount = 0;
    shadowMapsInfo.pQueueFamilyIndices = nullptr;
    shadowMapsInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

    m_shadowMaps = m_memoryManager.requestImage(m_renderDevice, shadowMapsInfo, VulkanMemoryType::TEXTURE);

    VkImageSubresourceRange subresourceRange;
    subresourceRange.aspectMask = VK_IMAGE_ASPECT_DEPTH_BIT;
    subresourceRange.baseMipLevel = 0;
    subresourceRange.levelCount = 1;
    subresourceRange.baseArrayLayer = 0;
    subresourceRange.layerCount = std::max(numShadows, (uint32_t)1);

    VkImageViewCreateInfo imageViewInfo;
    imageViewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
    imageViewInfo.pNext = nullptr;
    imageViewInfo.flags = 0;
    imageViewInfo.image = *m_shadowMaps->getImage();
    imageViewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D_ARRAY;
    imageViewInfo.format = VK_FORMAT_D32_SFLOAT;
    imageViewInfo.components = VulkanDefaults::getDefaultComponentMapping();
    imageViewInfo.subresourceRange = subresourceRange;

    vkCreateImageView(m_renderDevice, &imageViewInfo, nullptr, &m_shadowMapsView);

    m_shadowFramebuffers.clear();
    m_shadowMapsViews.resize((size_t)numShadows);
    m_shadowPasses.resize((size_t)numShadows);
    uint32_t currentLight = 0;
    auto shadowSamples = VK_SAMPLE_COUNT_1_BIT;

    for (auto light : m_scene->getLights())
    {
        if (light->getType() == LightType::DIRECTIONAL_LIGHT && currentLight < 16)
        {
            imageViewInfo.subresourceRange.baseArrayLayer = currentLight;
            imageViewInfo.subresourceRange.layerCount = 1;
            imageViewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;

            vkCreateImageView(m_renderDevice, &imageViewInfo, nullptr, &m_shadowMapsViews[currentLight]);

            VulkanRenderPassGenerator::generateShadowRenderPass(m_renderDevice, m_shadowPasses[currentLight], shadowSamples);

            m_shadowFramebuffers.push_back(
                std::make_shared<VulkanFramebuffer>(m_memoryManager, resolution, resolution, false, shadowSamples));
            m_shadowFramebuffers[currentLight]->setDepth(&m_shadowMapsViews[currentLight], VK_FORMAT_D32_SFLOAT);
            m_shadowFramebuffers[currentLight]->initializeFramebuffer(&m_shadowPasses[currentLight]);


            auto directionalLight = std::dynamic_pointer_cast<DirectionalLight>(light);
            directionalLight->m_shadowMapIndex = currentLight;
            m_shadowLights.push_back(directionalLight);
            currentLight++;
        }
    }

    m_lightMatrices.resize(currentLight);
}

void
VulkanRenderer::setShadowMapResolution(uint32_t resolution)
{
    m_shadowMapResolution = resolution;
}

void
VulkanRenderer::setResolution(unsigned int width, unsigned int height)
{
    m_width = width;
    m_height = height;
}

void
VulkanRenderer::setBloomOn()
{
    m_postProcessingChain->m_bloom = true;
}

void
VulkanRenderer::setBloomOff()
{
    m_postProcessingChain->m_bloom = false;
}

void
VulkanRenderer::setCommandBufferState(VkCommandBuffer * commandBuffer, uint32_t width, uint32_t height)
{
    VkViewport viewport;
    viewport.x = 0;
    viewport.y = 0;
    viewport.height = height;
    viewport.width = width;
    viewport.minDepth = 0.0;
    viewport.maxDepth = 1.0;

    vkCmdSetViewport(*commandBuffer, 0, 1, &viewport);

    VkRect2D scissor;

    scissor.offset = { 0, 0 };
    scissor.extent = { (uint32_t)viewport.width,
                       (uint32_t)viewport.height };

    vkCmdSetScissor(*commandBuffer, 0, 1, &scissor);
}

VulkanRenderer::~VulkanRenderer()
{
    // Delete devices
    for (int i = 0; i < (int)(m_deviceCount); i++)
    {
        // Important: must wait for device to be finished
        vkDeviceWaitIdle(m_devices[i]);

        vkDestroySemaphore(m_renderDevice, m_readyToRender, nullptr);
        vkDestroySemaphore(m_renderDevice, m_drawingComplete, nullptr);
        vkDestroySemaphore(m_renderDevice, m_presentImages, nullptr);

        for (auto fence : m_commandBufferSubmit)
        {
            vkDestroyFence(m_renderDevice, fence, nullptr);
        }

        // Clear all memory
        m_memoryManager.clear();

        // Delete framebuffers
        this->deleteFramebuffers();

        // Delete shadows
        for (auto imageView : m_shadowMapsViews)
        {
            vkDestroyImageView(m_renderDevice, imageView, nullptr);
        }
        vkDestroyImageView(m_renderDevice, m_shadowMapsView, nullptr);

        // Delete textures
        for (auto texture : m_textureMap)
        {
            texture.second->clear(&m_renderDevice);
        }
        m_noiseTextureDelegate->clear(&m_renderDevice);

        for (auto renderDelegate : m_renderDelegates)
        {
            renderDelegate->m_material->clear(&m_renderDevice);
            renderDelegate->m_depthMaterial->clear(&m_renderDevice);
            if (m_shadowPasses.size() > 0)
            {
                renderDelegate->m_shadowMaterial->clear(&m_renderDevice);
            }
        }

        for (auto postProcess : m_postProcessingChain->m_postProcesses)
        {
            postProcess->clear(&m_renderDevice);
        }

        for (auto pass : m_ssao)
        {
            pass->clear(&m_renderDevice);
        }

        for (auto pass : m_HDRTonemaps)
        {
            pass->clear(&m_renderDevice);
        }

        vkDestroyPipelineCache(m_renderDevice, m_pipelineCache, nullptr);

        vkDestroyRenderPass(m_renderDevice, m_opaqueRenderPass, nullptr);
        vkDestroyRenderPass(m_renderDevice, m_decalRenderPass, nullptr);
        vkDestroyRenderPass(m_renderDevice, m_depthRenderPass, nullptr);

        for (auto pass : m_shadowPasses)
        {
            vkDestroyRenderPass(m_renderDevice, pass, nullptr);
        }

        vkDestroySampler(m_renderDevice, m_HDRImageSampler, nullptr);

        for (auto sampler : m_swapchainImageSamplers)
        {
            vkDestroySampler(m_renderDevice, sampler, nullptr);
        }

        for (auto framebuffer : m_shadowFramebuffers)
        {
            framebuffer->clear(&m_renderDevice);
        }

        vkDestroySwapchainKHR(m_renderDevice, *m_swapchain, nullptr);

        // Delete command buffers
        vkDestroyCommandPool(m_renderDevice, m_renderCommandPool, nullptr);
        vkDestroyCommandPool(m_renderDevice, m_postProcessingCommandPool, nullptr);

        vkDestroyDevice(m_devices[i], nullptr);
    }

#ifndef NDEBUG
    PFN_vkDestroyDebugReportCallbackEXT createCallback =
        (PFN_vkDestroyDebugReportCallbackEXT)vkGetInstanceProcAddr(*m_instance, "vkDestroyDebugReportCallbackEXT");
    createCallback(*m_instance, m_debugReportCallback, nullptr);
#endif
    vkDestroyInstance(*m_instance, nullptr);
}
}