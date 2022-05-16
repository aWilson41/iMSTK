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

#pragma once

#include "imstkCellMesh.h"
#include "imstkMacros.h"

#include <array>
#include <unordered_set>

namespace imstk
{
struct Color;

///
/// \class LineMesh
///
/// \brief Base class for all volume mesh types
///
class LineMesh : public CellMesh<2>
{
public:
    LineMesh() = default;
    ~LineMesh() override = default;

    IMSTK_TYPE_NAME(LineMesh)

    ///
    /// \brief
    ///
    void clear() override;

    ///
    /// \brief
    ///
    void print() const override;

// Attributes
    ///
    /// \brief Get/Set the active scalars
    ///@{
    void setCellScalars(const std::string& arrayName, std::shared_ptr<AbstractDataArray> scalars);
    void setCellScalars(const std::string& arrayName);
    std::string getActiveCellScalars() const { return m_activeCellScalars; }
    std::shared_ptr<AbstractDataArray> getCellScalars() const;
///@}

protected:
    std::string m_activeCellScalars = "";
};
} // namespace imstk
