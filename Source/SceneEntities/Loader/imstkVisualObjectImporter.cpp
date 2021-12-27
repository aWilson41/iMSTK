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

#include "imstkVisualObjectImporter.h"
#include "imstkAnimatedObject.h"
#include "imstkAnimationModel.h"
#include "imstkAnimationNode.h"
#include "imstkAssimpMeshIO.h"
#include "imstkAttachment.h"
#include "imstkKeyFrameNodeAnimation.h"
#include "imstkLogger.h"
#include "imstkRenderMaterial.h"
#include "imstkSurfaceMesh.h"
#include "imstkVecDataArray.h"
#include "imstkVisualModel.h"

#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <assimp/scene.h>

#include <stack>

namespace imstk
{
static Mat4d
aiMatToMat4d(const aiMatrix4x4& m)
{
    Mat4d results;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            results(i, j) = static_cast<double>(m[i][j]);
        }
    }
    return results;
}

static aiMatrix4x4
mat4dToAiMat(const Mat4d& m)
{
    aiMatrix4x4 results;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            results[i][j] = static_cast<float>(m(i, j));
        }
    }
    return results;
}

std::shared_ptr<SceneObject>
ObjectIO::importSceneObject(
    const std::string& objName,
    const std::string& modelFilePath,
    const std::string& textureFolderPath,
    const Mat4d&       transform)
{
    auto type = MeshIO::getFileType(modelFilePath);

    // Check if mesh is supported by Assimp
    if (type != MeshFileType::_3DS
        && type != MeshFileType::OBJ
        && type != MeshFileType::FBX
        && type != MeshFileType::DAE)
    {
        LOG(FATAL) << "Error: File type not supported! Input model file path: " << modelFilePath;
        return nullptr;
    }

    auto visualObject = std::make_shared<SceneObject>(objName);

    // Import mesh(es) and apply some clean-up operations
    Assimp::Importer importer;
    const aiScene*   scene = importer.ReadFile(modelFilePath, AssimpMeshIO::getDefaultPostProcessSteps());

    // Check if there is actually a mesh or if the file can be read
    CHECK(scene != nullptr && scene->HasMeshes()) << "Error: could not read model with Assimp reader! Input model file path: "
                                                  << modelFilePath;

    // Iterate over each material
    std::vector<std::shared_ptr<RenderMaterial>> materials(scene->mNumMaterials);
    for (unsigned int i = 0; i < scene->mNumMaterials; i++)
    {
        materials[i] = readMaterial(scene->mMaterials[i], textureFolderPath);
    }

    // Read all meshes
    std::vector<std::shared_ptr<PointSet>>       meshes(scene->mNumMeshes);
    std::vector<std::shared_ptr<RenderMaterial>> meshMaterials(scene->mNumMeshes);
    for (unsigned int i = 0; i < scene->mNumMeshes; i++)
    {
        aiMesh* importedMesh = scene->mMeshes[i];
        meshes[i] = AssimpMeshIO::convertAssimpMesh(importedMesh);
        meshMaterials[i] = materials[importedMesh->mMaterialIndex];
    }

    // Iterate through assimp graph
    std::stack<aiNode*> nodes;
    nodes.push(scene->mRootNode);
    scene->mRootNode->mTransformation = mat4dToAiMat(transform) * scene->mRootNode->mTransformation;
    std::unordered_map<aiNode*, Mat4d> worldTransforms;
    while (!nodes.empty())
    {
        aiNode* currNode = nodes.top();
        nodes.pop();

        Mat4d parentWorldTransform = Mat4d::Identity();
        if (currNode->mParent != nullptr)
        {
            parentWorldTransform = worldTransforms[currNode->mParent];
        }
        Mat4d localTransform     = aiMatToMat4d(currNode->mTransformation);
        Mat4d currWorldTransform = parentWorldTransform * localTransform;
        worldTransforms[currNode] = currWorldTransform;

        for (unsigned int i = 0; i < currNode->mNumMeshes; i++)
        {
            // Copy, transform, and insert the mesh
            std::shared_ptr<SurfaceMesh> surfMesh  = std::make_shared<SurfaceMesh>();
            const unsigned int           meshIndex = currNode->mMeshes[i];
            surfMesh->deepCopy(std::dynamic_pointer_cast<SurfaceMesh>(meshes[meshIndex]));
            auto visualModel = std::make_shared<VisualModel>(surfMesh);
            visualModel->setName(std::string(currNode->mName.C_Str()));

            surfMesh->transform(currWorldTransform, Geometry::TransformType::ApplyToData);
            visualModel->setRenderMaterial(meshMaterials[meshIndex]);
            visualObject->addVisualModel(visualModel);
        }

        for (unsigned int i = 0; i < currNode->mNumChildren; i++)
        {
            nodes.push(currNode->mChildren[i]);
        }
    }

    return visualObject;
}

std::shared_ptr<AnimatedObject>
ObjectIO::importAnimObject(
    const std::string& objName,
    const std::string& modelFilePath,
    const std::string& textureFolderPath,
    const Mat4d&       transform)
{
    const MeshFileType type = MeshIO::getFileType(modelFilePath);
    if (type != MeshFileType::_3DS
        && type != MeshFileType::FBX
        && type != MeshFileType::DAE
        && type != MeshFileType::GLTF
        && type != MeshFileType::BLEND
        && type != MeshFileType::X)
    {
        LOG(FATAL) << "File type not supported";
        return nullptr;
    }

    // Import the scene
    Assimp::Importer importer;
    const aiScene*   scene = importer.ReadFile(modelFilePath,
        aiPostProcessSteps::aiProcess_GenSmoothNormals |
        aiPostProcessSteps::aiProcess_CalcTangentSpace |
        aiPostProcessSteps::aiProcess_Triangulate |
        aiPostProcessSteps::aiProcess_ImproveCacheLocality);
    ;
    CHECK(scene != nullptr) << "AssimpMeshIO::readMeshData error: Failed to read "
        << objName << " from " << modelFilePath << " " << importer.GetErrorString();

    // Read every material (they could be shared later)
    std::vector<std::shared_ptr<RenderMaterial>> materials(scene->mNumMaterials);
    for (unsigned int i = 0; i < scene->mNumMaterials; i++)
    {
        materials[i] = readMaterial(scene->mMaterials[i], textureFolderPath);
    }

    // Read the nodes, establish some maps for later
    auto rootNode = std::make_shared<AnimationNode>();

    std::unordered_map<std::string, std::shared_ptr<AnimationNode>> nameToNodeMap;
    std::unordered_map<std::shared_ptr<AnimationNode>, aiNode*>     nodeToAiNodeMap;
    std::unordered_map<aiNode*, std::shared_ptr<AnimationNode>>     aiNodeToNodeMap;
    std::unordered_map<int, aiNode*>                                meshIndexToAiNodeMap;
    {
        // Setup a stack to DFS down the tree
        std::stack<aiNode*> aiNodes;
        aiNodes.push(scene->mRootNode);

        aiNodeToNodeMap[scene->mRootNode] = rootNode;

        while (!aiNodes.empty())
        {
            aiNode* aiCurrNode = aiNodes.top();
            aiNodes.pop();

            // Get current node from the map
            std::shared_ptr<AnimationNode> currNode = aiNodeToNodeMap[aiCurrNode];
            currNode->name = std::string(aiCurrNode->mName.C_Str());
            currNode->localBindTransform  = aiMatToMat4d(aiCurrNode->mTransformation);
            nameToNodeMap[currNode->name] = currNode;
            nodeToAiNodeMap[currNode]     = aiCurrNode;

            for (unsigned int i = 0; i < aiCurrNode->mNumMeshes; i++)
            {
                meshIndexToAiNodeMap[aiCurrNode->mMeshes[i]] = aiCurrNode;
            }

            for (unsigned int i = 0; i < aiCurrNode->mNumChildren; i++)
            {
                // Get the assimp child
                aiNode* aiChildNode = aiCurrNode->mChildren[i];

                // Create and link the child
                auto childNode = std::make_shared<AnimationNode>();
                currNode->children.push_back(childNode);
                childNode->parent = currNode.get();

                // Map the child
                aiNodeToNodeMap[aiChildNode] = childNode;

                aiNodes.push(aiCurrNode->mChildren[i]);
            }
        }

        // Edge case, we have an additional transform for the root
        rootNode->localBindTransform = transform * aiMatToMat4d(scene->mRootNode->mTransformation);
    }

    std::unordered_map<std::shared_ptr<Geometry>, std::unordered_map<int, double>> vertexWeightSums; // For later normalization

    // Read all meshes (setup bone weights as well)
    std::vector<std::shared_ptr<SurfaceMesh>>    meshes(scene->mNumMeshes);
    std::vector<std::shared_ptr<RenderMaterial>> meshMaterials(scene->mNumMeshes);
    for (unsigned int i = 0; i < scene->mNumMeshes; i++)
    {
        aiMesh* importedMesh = scene->mMeshes[i];
        meshes[i] = std::dynamic_pointer_cast<SurfaceMesh>(AssimpMeshIO::convertAssimpMesh(importedMesh));
        meshMaterials[i] = materials[importedMesh->mMaterialIndex];

        // We need to store pre-post normals for skeletal animations
        // The normal needs to be rotated in the shader (or in our case on the CPU)
        std::shared_ptr<VecDataArray<double, 3>> normals = meshes[i]->getVertexNormals();
        if (normals != nullptr)
        {
            meshes[i]->setVertexAttribute("preNormals", std::make_shared<VecDataArray<double, 3>>(*normals));
        }

        // Identify all AnimationNodes that modify this mesh
        for (unsigned int j = 0; j < importedMesh->mNumBones; j++)
        {
            aiBone*                        bone     = importedMesh->mBones[j];
            std::shared_ptr<AnimationNode> boneNode = nameToNodeMap[std::string(bone->mName.C_Str())];
            boneNode->boneOffset = aiMatToMat4d(bone->mOffsetMatrix);

            if (bone->mNumWeights > 0)
            {
                // If the bone has weights, then add the mesh as a MeshWeightAttachment
                // to do weighted transforms with the node
                auto meshWeightAttachment = std::make_shared<MeshWeightAttachment>();
                meshWeightAttachment->geometry = meshes[i];
                meshWeightAttachment->weights.resize(bone->mNumWeights);

                for (unsigned int k = 0; k < bone->mNumWeights; k++)
                {
                    VertexWeight weight;
                    weight.vertexId = bone->mWeights[k].mVertexId;
                    weight.weight   = bone->mWeights[k].mWeight;

                    meshWeightAttachment->weights[k] = weight;

                    vertexWeightSums[meshWeightAttachment->geometry][weight.vertexId];
                    vertexWeightSums[meshWeightAttachment->geometry][weight.vertexId] += weight.weight;
                }

                meshWeightAttachment->parent = boneNode.get();
                boneNode->attachments.push_back(meshWeightAttachment);
            }
        }
        // If it has no bones add the mesh as a standard attachment
        // This will transform it with the node
        if (importedMesh->mNumBones == 0)
        {
            std::shared_ptr<AnimationNode> node = aiNodeToNodeMap[meshIndexToAiNodeMap[i]];
            auto                           meshAttachment = std::make_shared<MeshAttachment>();
            meshAttachment->geometry = meshes[i];
            meshAttachment->parent   = node.get();
            node->attachments.push_back(meshAttachment);
        }
    }

    // Normalize the weights (lot of file types don't have this so when applying a transform that
    // isn't what's given in the file we get weird squiggly transforms)
    for (auto i : nameToNodeMap)
    {
        std::shared_ptr<AnimationNode> node = i.second;
        for (auto attachment : node->attachments)
        {
            if (auto weightAttachment = std::dynamic_pointer_cast<MeshWeightAttachment>(attachment))
            {
                for (VertexWeight& weight : weightAttachment->weights)
                {
                    weight.weight /= vertexWeightSums[weightAttachment->geometry][weight.vertexId];
                }
            }
        }
    }

    auto animModel = std::make_shared<AnimationModel>();
    animModel->setRootNode(rootNode);

    // Read the animations
    for (unsigned int i = 0; i < scene->mNumAnimations; i++)
    {
        aiAnimation*                    anim      = scene->mAnimations[i];
        std::shared_ptr<AnimationTrack> animation = std::make_shared<AnimationTrack>();
        animation->channels.resize(anim->mNumChannels);
        //printf("Num channels %d\n", anim->mNumChannels);
        for (unsigned int j = 0; j < anim->mNumChannels; j++)
        {
            aiNodeAnim*                    aiNodeAnim = anim->mChannels[j];
            const std::string nodeName = std::string(aiNodeAnim->mNodeName.C_Str());
            std::shared_ptr<AnimationNode> node = nameToNodeMap[nodeName];

            std::shared_ptr<KeyFrameNodeAnimation> nodeAnim   = std::make_shared<KeyFrameNodeAnimation>();
            nodeAnim->setNode(node);
            //printf("Num Keys: %d\n", aiNodeAnim->mNumPositionKeys);
            for (unsigned int k = 0; k < aiNodeAnim->mNumPositionKeys; k++)
            {
                aiVectorKey keyValue = aiNodeAnim->mPositionKeys[k];
                nodeAnim->m_positions[keyValue.mTime] =
                    Vec3d(keyValue.mValue.x, keyValue.mValue.y, keyValue.mValue.z);
            }
            for (unsigned int k = 0; k < aiNodeAnim->mNumRotationKeys; k++)
            {
                aiQuatKey keyValue = aiNodeAnim->mRotationKeys[k];
                nodeAnim->m_orientations[keyValue.mTime] =
                    Quatd(keyValue.mValue.w, keyValue.mValue.x, keyValue.mValue.y, keyValue.mValue.z);
            }
            for (unsigned int k = 0; k < aiNodeAnim->mNumScalingKeys; k++)
            {
                aiVectorKey keyValue = aiNodeAnim->mScalingKeys[k];
                nodeAnim->m_scalings[keyValue.mTime] =
                    Vec3d(keyValue.mValue.x, keyValue.mValue.y, keyValue.mValue.z);
            }
            animation->channels[j] = nodeAnim;
        }
        animModel->addAnimation(std::string(anim->mName.C_Str()), animation);
        //printf("Animation %s added\n", anim->mName.C_Str());
    }

    // Finally create the object and start setting it up
    animModel->setTime(0.0);

    auto animObject = std::make_shared<AnimatedObject>(objName);
    animObject->setAnimationModel(animModel);

    // Note: We do not support instancing. So if a mesh is referenced twice in the tree,
    // it will fail. Every mesh is simply added once
    // or else we'd need to import multiple AnimatedObject's
    for (int i = 0; i < meshes.size(); i++)
    {
        auto visualModel = std::make_shared<VisualModel>(meshes[i]);
        visualModel->setRenderMaterial(meshMaterials[i]);
        visualModel->getRenderMaterial()->setRecomputeVertexNormals(false);
        animObject->addVisualModel(visualModel);
    }

    return animObject;
}

std::shared_ptr<Texture>
ObjectIO::createTexture(std::string textureFolderPath, std::string textureFilePath, Texture::Type textureType)
{
    textureFilePath = getSubstringGivenString(textureFilePath, "/", true);
    textureFilePath = getSubstringGivenString(textureFilePath, "\\", true);

    std::string fileName = getSubstringGivenString(textureFilePath, ".", false);

    const std::string fileExt = getSubstringGivenString(textureFilePath, ".", true);

    const std::string filePath = textureFolderPath + fileName + "." + fileExt;

    // Check if file exists
    std::ifstream file(filePath);
    if (file.good())
    {
        file.close();
        return std::make_shared<Texture>(filePath, textureType);
    }
    return nullptr;
}

std::shared_ptr<RenderMaterial>
ObjectIO::readMaterial(aiMaterial* material, std::string textureFolderPath)
{
    // todo: Embedded textures (GLTF2.0)

    // Create our material
    auto renderMaterial = std::make_shared<RenderMaterial>();
    renderMaterial->setShadingModel(RenderMaterial::ShadingModel::Phong);

    // Go through all the ai properties, not all included here
    aiString name;
    aiReturn ret = material->Get(AI_MATKEY_NAME, name);
    if (ret == AI_SUCCESS)
    {
        renderMaterial->setName(std::string(name.C_Str()));
    }

    aiColor3D ambientColor;
    ret = material->Get(AI_MATKEY_COLOR_AMBIENT, ambientColor);
    if (ret == AI_SUCCESS)
    {
        renderMaterial->setAmbientColor(Color(ambientColor.r, ambientColor.g, ambientColor.b));
    }

    aiColor3D diffuseColor;
    ret = material->Get(AI_MATKEY_COLOR_DIFFUSE, diffuseColor);
    if (ret == AI_SUCCESS)
    {
        renderMaterial->setDiffuseColor(Color(diffuseColor.r, diffuseColor.g, diffuseColor.b));
    }

    aiColor3D specularColor;
    ret = material->Get(AI_MATKEY_COLOR_SPECULAR, specularColor);
    if (ret == AI_SUCCESS)
    {
        renderMaterial->setSpecularColor(Color(specularColor.r, specularColor.g, specularColor.b));
    }

    int useWireframe;
    ret = material->Get(AI_MATKEY_ENABLE_WIREFRAME, useWireframe);
    if (ret == AI_SUCCESS)
    {
        if (useWireframe == 1)
        {
            renderMaterial->setDisplayMode(RenderMaterial::DisplayMode::Wireframe);
        }
    }

    int useBackfaceRendering;
    ret = material->Get(AI_MATKEY_TWOSIDED, useBackfaceRendering);
    if (ret == AI_SUCCESS)
    {
        renderMaterial->setBackFaceCulling(static_cast<bool>(useBackfaceRendering));
    }

    float opacity;
    ret = material->Get(AI_MATKEY_OPACITY, opacity);
    if (ret == AI_SUCCESS)
    {
        renderMaterial->setOpacity(opacity);
    }

    float shininess;
    ret = material->Get(AI_MATKEY_SHININESS, shininess);
    if (ret == AI_SUCCESS)
    {
        renderMaterial->setSpecular(shininess);
    }

    float shininessStrength;
    ret = material->Get(AI_MATKEY_SHININESS_STRENGTH, shininessStrength);
    if (ret == AI_SUCCESS)
    {
        renderMaterial->setSpecularPower(shininessStrength);
    }

    float reflectivity;
    ret = material->Get(AI_MATKEY_REFLECTIVITY, reflectivity);
    if (ret == AI_SUCCESS)
    {
        // Not supported yet
    }

    renderMaterial->setRecomputeVertexNormals(false);

    aiString texFilePath;
    ret = material->GetTexture(aiTextureType::aiTextureType_AMBIENT, 0, &texFilePath);
    if (ret == AI_SUCCESS)
    {
        std::shared_ptr<Texture> tex = createTexture(textureFolderPath, std::string(texFilePath.C_Str()), Texture::Type::AmbientOcclusion);
        if (tex != nullptr)
        {
            renderMaterial->addTexture(tex);
        }
    }
    ret = material->GetTexture(aiTextureType::aiTextureType_DIFFUSE, 0, &texFilePath);
    if (ret == AI_SUCCESS)
    {
        std::shared_ptr<Texture> tex = createTexture(textureFolderPath, std::string(texFilePath.C_Str()), Texture::Type::Diffuse);
        if (tex != nullptr)
        {
            renderMaterial->addTexture(tex);
        }
    }
    ret = material->GetTexture(aiTextureType::aiTextureType_EMISSIVE, 0, &texFilePath);
    if (ret == AI_SUCCESS)
    {
        std::shared_ptr<Texture> tex = createTexture(textureFolderPath, std::string(texFilePath.C_Str()), Texture::Type::Emissive);
        if (tex != nullptr)
        {
            renderMaterial->addTexture(tex);
        }
    }
    ret = material->GetTexture(aiTextureType::aiTextureType_NORMALS, 0, &texFilePath);
    if (ret == AI_SUCCESS)
    {
        std::shared_ptr<Texture> tex = createTexture(textureFolderPath, std::string(texFilePath.C_Str()), Texture::Type::Normal);
        if (tex != nullptr)
        {
            renderMaterial->addTexture(tex);
        }
    }
    ret = material->GetTexture(aiTextureType::aiTextureType_SPECULAR, 0, &texFilePath);
    if (ret == AI_SUCCESS)
    {
        std::shared_ptr<Texture> tex = createTexture(textureFolderPath, std::string(texFilePath.C_Str()), Texture::Type::Metalness);
        if (tex != nullptr)
        {
            renderMaterial->addTexture(tex);
        }
    }
    return renderMaterial;
}

std::string
ObjectIO::getSubstringGivenString(
    const std::string& input,
    const std::string& delimiter,
    const bool         lastInstance /*= false*/)
{
    unsigned long long index = 0;
    unsigned long long tempIndex;

    if (lastInstance)
    {
        tempIndex = input.rfind(delimiter) + 1;
        if (tempIndex >= input.length())
        {
            return input;
        }
    }
    else
    {
        tempIndex = input.find(delimiter);
    }

    if (tempIndex == std::string::npos)
    {
        index = lastInstance ? 0 : input.length();
    }
    else
    {
        index = tempIndex;
    }

    if (lastInstance)
    {
        return input.substr(index);
    }

    return input.substr(0, index);
}
} // imstk
