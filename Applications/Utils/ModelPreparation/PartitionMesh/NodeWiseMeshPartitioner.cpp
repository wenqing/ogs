/*!
  \file NodeWiseMeshPartitioner.cpp
  \date   2016.05

  \brief  Define the members of class NodeWiseMeshPartitioner

  \copyright
  Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include "NodeWiseMeshPartitioner.h"

#include <limits>
#include <iomanip>
#include <cstdio>  // for binary output
#include <numeric>

#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"

#include "MeshLib/IO/VtkIO/VtuInterface.h"

#include "MeshLib/Elements/Element.h"

namespace ApplicationUtils
{
struct NodeStruct
{
    NodeWiseMeshPartitioner::IntegerType id;
    double x;
    double y;
    double z;
};

void NodeWiseMeshPartitioner::readMetisData(const std::string& file_name_base)
{
    const std::string npartitions_str = std::to_string(_npartitions);

    // Read partitioned mesh data from METIS
    const std::string fname_parts = file_name_base + ".mesh.npart." + npartitions_str;

    std::ifstream npart_in(fname_parts);
    if (!npart_in.is_open())
    {
        OGS_FATAL(
            "Error: cannot open file %s. It may not exist!\n"
            "Run mpmetis beforehand or use option -m",
            fname_parts.data());
    }

    const std::size_t nnodes = _mesh->getNumberOfNodes();

    std::size_t counter = 0;
    while (!npart_in.eof())
    {
        npart_in >> _nodes_partition_ids[counter++] >> std::ws;
        if (counter == nnodes)
          break;
    }

    if (npart_in.bad())
    {
        OGS_FATAL(
            "Error while reading file %s.", fname_parts.data());
    }

    npart_in.close();

    if( counter != nnodes)
    {
        OGS_FATAL(
            "Error: data in %s are less than expected.", fname_parts.data());
    }

    // remove metis files.
    std::remove(fname_parts.c_str());
    const std::string fname_eparts = file_name_base + ".mesh.epart."
                                     + npartitions_str;
    std::remove(fname_eparts.c_str());
}

void NodeWiseMeshPartitioner::findNonGhostNodesInPartition(
    std::size_t const part_id,
    const bool is_mixed_high_order_linear_elems,
    std::vector<MeshLib::Node*>& extra_nodes)
{
    std::vector<MeshLib::Node*> const& nodes = _mesh->getNodes();
    auto& partition = _partitions[part_id];
    // -- Extra nodes for high order elements
    for (std::size_t i = 0; i < _mesh->getNumberOfNodes(); i++)
    {
        if (_nodes_partition_ids[i] == part_id)
        {
            splitOfHigherOrderNode(nodes, is_mixed_high_order_linear_elems,
                                   i, partition.nodes, extra_nodes);
        }
    }
    partition.number_of_non_ghost_base_nodes = partition.nodes.size();
    partition.number_of_non_ghost_nodes =
        partition.number_of_non_ghost_base_nodes + extra_nodes.size();
}

void NodeWiseMeshPartitioner::findElementsInPartition(
    std::size_t const part_id, const bool is_mixed_high_order_linear_elems)
{
    auto& partition = _partitions[part_id];
    std::vector<MeshLib::Element*> const& elements = _mesh->getElements();
    std::vector<bool> _is_regular_element(elements.size(), false);

    for (std::size_t elem_id = 0; elem_id < elements.size(); elem_id++)
    {
        const auto* elem = elements[elem_id];
        if (_is_regular_element[elem_id])
            continue;

        std::size_t non_ghost_node_number = 0;
        for (unsigned i = 0; i < elem->getNumberOfNodes(); i++)
        {
            if (_nodes_partition_ids[elem->getNodeIndex(i)] == part_id)
            {
                non_ghost_node_number++;
            }
        }

        if (non_ghost_node_number == 0)
            continue;

        if (non_ghost_node_number == elem->getNumberOfNodes())
        {
            partition.regular_elements.push_back(elem);
            _is_regular_element[elem_id] = true;
        }
        else
        {
            partition.ghost_elements.push_back(elem);
        }
    }
}

void NodeWiseMeshPartitioner::findGhostNodesInPartition(
    std::size_t const part_id,
    const bool is_mixed_high_order_linear_elems,
    std::vector<MeshLib::Node*>& extra_nodes)
{
    auto& partition = _partitions[part_id];
    std::vector<MeshLib::Node*> const& nodes = _mesh->getNodes();
    std::vector<bool> nodes_reserved(_mesh->getNumberOfNodes(), false);
    for (const auto* ghost_elem : partition.ghost_elements)
    {
        for (unsigned i = 0; i < ghost_elem->getNumberOfNodes(); i++)
        {
            const unsigned node_id = ghost_elem->getNodeIndex(i);
            if (nodes_reserved[node_id])
                continue;

            if (_nodes_partition_ids[node_id] != part_id)
            {
                splitOfHigherOrderNode(nodes, is_mixed_high_order_linear_elems,
                                       node_id, partition.nodes, extra_nodes);
                nodes_reserved[node_id] = true;
            }
        }
    }
}

void NodeWiseMeshPartitioner::splitOfHigherOrderNode(
    std::vector<MeshLib::Node*> const& nodes,
    bool const is_mixed_high_order_linear_elems,
    unsigned const node_id,
    std::vector<MeshLib::Node*>& base_nodes,
    std::vector<MeshLib::Node*>& extra_nodes)
{
    if (is_mixed_high_order_linear_elems)
    {
        if (node_id < _mesh->getNumberOfBaseNodes())
            base_nodes.push_back(nodes[node_id]);
        else
            extra_nodes.push_back(nodes[node_id]);
    }
    else
    {
        base_nodes.push_back(nodes[node_id]);
    }
}

void NodeWiseMeshPartitioner::processPartition(std::size_t const part_id,
    const bool is_mixed_high_order_linear_elems)
{
    std::vector<MeshLib::Node*> extra_nodes;
    findNonGhostNodesInPartition(part_id, is_mixed_high_order_linear_elems,
                                 extra_nodes);

    findElementsInPartition(part_id, is_mixed_high_order_linear_elems);
    findGhostNodesInPartition(part_id, is_mixed_high_order_linear_elems,
                              extra_nodes);
    auto& partition = _partitions[part_id];
    partition.number_of_base_nodes = partition.nodes.size();

    if (is_mixed_high_order_linear_elems)
        partition.nodes.insert(partition.nodes.end(), extra_nodes.begin(),
                               extra_nodes.end());
}

void NodeWiseMeshPartitioner::processProperties()
{
    std::size_t const total_number_of_tuples =
        std::accumulate(std::begin(_partitions), std::end(_partitions), 0,
                        [](std::size_t const sum, Partition const& p) {
                            return sum + p.nodes.size();
                        });

    INFO("total number of tuples after partitioning: %d ",
         total_number_of_tuples);
    // 1 create new PV
    // 2 resize the PV with total_number_of_tuples
    // 3 copy the values according to the partition info
    auto const& original_properties(_mesh->getProperties());
    auto property_names = original_properties.getPropertyVectorNames();
    for (auto const& name : property_names)
    {
        bool success =
            copyPropertyVector<double>(name, total_number_of_tuples) ||
            copyPropertyVector<float>(name, total_number_of_tuples) ||
            copyPropertyVector<int>(name, total_number_of_tuples) ||
            copyPropertyVector<long>(name, total_number_of_tuples) ||
            copyPropertyVector<unsigned>(name, total_number_of_tuples) ||
            copyPropertyVector<unsigned long>(name, total_number_of_tuples) ||
            copyPropertyVector<std::size_t>(name, total_number_of_tuples);
        if (!success)
            WARN(
                "processProperties: Could not create partitioned "
                "PropertyVector '%s'.",
                name.c_str());
    }
}

void NodeWiseMeshPartitioner::partitionByMETIS(
    const bool is_mixed_high_order_linear_elems)
{
    for (std::size_t part_id = 0; part_id < _partitions.size(); part_id++)
    {
        INFO("Processing partition: %d", part_id);
        processPartition(part_id, is_mixed_high_order_linear_elems);
    }

    renumberNodeIndices(is_mixed_high_order_linear_elems);

    processProperties();
}

void NodeWiseMeshPartitioner::renumberNodeIndices(
    const bool is_mixed_high_order_linear_elems)
{
    std::size_t node_global_id_offset = 0;
    // Renumber the global indices.
    // -- Base nodes
    for (auto& partition : _partitions)
    {
        for (std::size_t i = 0; i < partition.number_of_non_ghost_base_nodes;
             i++)
        {
            _nodes_global_ids[partition.nodes[i]->getID()] =
                node_global_id_offset;
            node_global_id_offset++;
        }
    }

    if (!is_mixed_high_order_linear_elems)
        return;

    // -- Nodes for high order elements.
    for (auto& partition : _partitions)
    {
        const std::size_t end_id = partition.number_of_base_nodes +
                                   partition.number_of_non_ghost_nodes -
                                   partition.number_of_non_ghost_base_nodes;
        for (std::size_t i = partition.number_of_base_nodes; i < end_id; i++)
        {
            _nodes_global_ids[partition.nodes[i]->getID()] =
                node_global_id_offset;
            node_global_id_offset++;
        }
    }
}

void NodeWiseMeshPartitioner::writeMETIS(const std::string& file_name)
{
    std::ofstream os(file_name, std::ios::trunc);
    if (!os.is_open())
    {
        OGS_FATAL("Error: cannot open file %s.",
                  file_name.data());
    }

    if (!os.good())
    {
        OGS_FATAL("Error: Cannot write in file %s.", file_name.data());
    }

    std::vector<MeshLib::Element*> const& elements = _mesh->getElements();
    os << elements.size() << " \n";
    for (const auto* elem : elements)
    {
        os << elem->getNodeIndex(0) + 1;
        for (unsigned j = 1; j < elem->getNumberOfNodes(); j++)
        {
            os << " " << elem->getNodeIndex(j) + 1;
        }
        os << "\n";
    }
}

NodeWiseMeshPartitioner::IntegerType
NodeWiseMeshPartitioner::getNumberOfIntegerVariablesOfElements(
    const std::vector<const MeshLib::Element*>& elements) const
{
    // Element ID, element type, and number of the nodes of
    // an element of all elements in the current partition.
    IntegerType nmb_element_idxs = 3 * elements.size();
    for (const auto* elem : elements)
    {
        nmb_element_idxs += elem->getNumberOfNodes();
    }
    return nmb_element_idxs;
}

void NodeWiseMeshPartitioner::writePropertiesBinary(
    const std::string& file_name_base) const
{
    auto const& property_names(_partitioned_properties.getPropertyVectorNames());
    if (property_names.empty())
        return;
    const std::string fname_cfg = file_name_base +
                                  "_partitioned_properties_cfg" +
                                  std::to_string(_npartitions) + ".bin";
    std::ofstream out(fname_cfg.c_str(), std::ios::binary | std::ios::out);

    const std::string fname_val = file_name_base +
                                  "_partitioned_properties_val" +
                                  std::to_string(_npartitions) + ".bin";
    std::ofstream out_val(fname_val.c_str(), std::ios::binary | std::ios::out);

    std::size_t number_of_properties(property_names.size());
    out.write(reinterpret_cast<char*>(&number_of_properties),
              sizeof(number_of_properties));
    for (auto const& name : property_names)
    {
        bool success =
            writePropertyVectorBinary<double>(name, out_val, out) ||
            writePropertyVectorBinary<float>(name, out_val, out) ||
            writePropertyVectorBinary<int>(name, out_val, out) ||
            writePropertyVectorBinary<long>(name, out_val, out) ||
            writePropertyVectorBinary<unsigned>(name, out_val, out) ||
            writePropertyVectorBinary<unsigned long>(name, out_val, out) ||
            writePropertyVectorBinary<std::size_t>(name, out_val, out);
        if (!success)
            OGS_FATAL(
                "writePropertiesBinary: Could not write PropertyVector '%s'.",
                name.c_str());
    }
    unsigned long offset = 0;
    for (const auto& partition : _partitions)
    {
        MeshLib::IO::PropertyVectorPartitionMetaData pvpmd;
        pvpmd.offset = offset;
        pvpmd.number_of_tuples = partition.nodes.size();
        INFO(
            "Write meta data for PropertyVector: global offset %d, number "
            "of tuples %d",
            pvpmd.offset, pvpmd.number_of_tuples);
        MeshLib::IO::writePropertyVectorPartitionMetaData(out, pvpmd);
        offset += pvpmd.number_of_tuples;
    }
    out.close();
    out_val.close();
}

std::tuple<std::vector<NodeWiseMeshPartitioner::IntegerType>,
           std::vector<NodeWiseMeshPartitioner::IntegerType>>
NodeWiseMeshPartitioner::writeConfigDataBinary(
    const std::string& file_name_base)
{
    const std::string fname = file_name_base + "_partitioned_msh_cfg"
                              + std::to_string(_npartitions) + ".bin";
    FILE* of_bin_cfg = fopen(fname.c_str(), "wb");

    const IntegerType num_config_data = 14;
    IntegerType config_data[num_config_data];
    // node rank offset
    config_data[10] = 0;
    // element rank offset
    config_data[11] = 0;
    // ghost element rank offset
    config_data[12] = 0;
    // Reserved
    config_data[13] = 0;
    std::vector<IntegerType> num_elem_integers(_partitions.size());
    std::vector<IntegerType> num_g_elem_integers(_partitions.size());
    std::size_t loop_id = 0;
    for (const auto& partition : _partitions)
    {
        config_data[0] = partition.nodes.size();
        config_data[1] = partition.number_of_base_nodes;
        config_data[2] = partition.regular_elements.size();
        config_data[3] = partition.ghost_elements.size();
        config_data[4] = partition.number_of_non_ghost_base_nodes;
        config_data[5] = partition.number_of_non_ghost_nodes;
        config_data[6] = _mesh->getNumberOfBaseNodes();
        config_data[7] = _mesh->getNumberOfNodes();
        config_data[8] =
            getNumberOfIntegerVariablesOfElements(partition.regular_elements);
        config_data[9] =
            getNumberOfIntegerVariablesOfElements(partition.ghost_elements);

        fwrite(config_data, 1, num_config_data * sizeof(IntegerType), of_bin_cfg);

        config_data[10] += config_data[0] * sizeof(NodeStruct);

        // Update offsets
        num_elem_integers[loop_id] =
            partition.regular_elements.size() + config_data[8];
        // Offset the ending entry of the element integer variales of
        // the non-ghost elements of this partition in the vector of elem_info.
        config_data[11] += num_elem_integers[loop_id] * sizeof(IntegerType);
        // Offset the ending entry of the element integer variales of
        // the ghost elements of this partition in the vector of elem_info.
        num_g_elem_integers[loop_id] =
            partition.ghost_elements.size() + config_data[9];
        config_data[12] += num_g_elem_integers[loop_id] * sizeof(IntegerType);

        loop_id++;
    }

    fclose(of_bin_cfg);

    return  std::make_tuple(num_elem_integers, num_g_elem_integers);
}

void NodeWiseMeshPartitioner::writeElementsBinary
                     (const std::string& file_name_base,
                      const std::vector<IntegerType>& num_elem_integers,
                      const std::vector<IntegerType>& num_g_elem_integers)
{
    const std::string npartitions_str = std::to_string(_npartitions);
    std::string fname = file_name_base + "_partitioned_msh_ele"
                              + npartitions_str + ".bin";
    FILE* of_bin_ele = fopen(fname.c_str(), "wb");
    fname =
        file_name_base + "_partitioned_msh_ele_g" + npartitions_str + ".bin";
    FILE* of_bin_ele_g = fopen(fname.c_str(), "wb");
    for (std::size_t i = 0; i < _partitions.size(); i++)
    {
        const auto& partition = _partitions[i];

        // Set the local node indices of the current partition.
        IntegerType node_local_id_offset = 0;
        std::vector<IntegerType> nodes_local_ids(_mesh->getNumberOfNodes(), -1);
        for (const auto* node : partition.nodes)
        {
            nodes_local_ids[node->getID()] = node_local_id_offset;
            node_local_id_offset++;
        }

        // A vector contians all element integer variales of
        // the non-ghost elements of this partition
        std::vector<IntegerType> ele_info(num_elem_integers[i]);

        // Non-ghost elements.
        IntegerType counter = partition.regular_elements.size();

        for (std::size_t j = 0; j < partition.regular_elements.size(); j++)
        {
            const auto* elem = partition.regular_elements[j];
            ele_info[j] = counter;
            getElementIntegerVariables(*elem, nodes_local_ids, ele_info,
                                       counter);
        }
        // Write vector data of non-ghost elements
        fwrite(ele_info.data(), 1, (num_elem_integers[i]) * sizeof(IntegerType),
               of_bin_ele);

        // Ghost elements
        ele_info.resize(num_g_elem_integers[i]);

        counter = partition.ghost_elements.size();

        for (std::size_t j = 0; j < partition.ghost_elements.size(); j++)
        {
            const auto* elem = partition.ghost_elements[j];
            ele_info[j] = counter;
            getElementIntegerVariables(*elem, nodes_local_ids, ele_info,
                                       counter);
        }
        // Write vector data of ghost elements
        fwrite(ele_info.data(), 1, (num_g_elem_integers[i]) * sizeof(IntegerType),
               of_bin_ele_g);
    }

    fclose(of_bin_ele);
    fclose(of_bin_ele_g);
}

void NodeWiseMeshPartitioner::writeNodesBinary(const std::string& file_name_base)
{
    const std::string fname = file_name_base + "_partitioned_msh_nod"
                              + std::to_string(_npartitions) + ".bin";
    FILE* of_bin_nod = fopen(fname.c_str(), "wb");

    for (const auto& partition : _partitions)
    {
        std::vector<NodeStruct> nodes_buffer;
        nodes_buffer.reserve(partition.nodes.size());

        for (const auto* node : partition.nodes)
        {
            double const* coords = node->getCoords();
            NodeStruct node_struct;
            node_struct.id = _nodes_global_ids[node->getID()];
            node_struct.x = coords[0];
            node_struct.y = coords[1];
            node_struct.z = coords[2];
            nodes_buffer.emplace_back(node_struct);
        }
        fwrite(nodes_buffer.data(), sizeof(NodeStruct), partition.nodes.size(),
               of_bin_nod);
    }
    fclose(of_bin_nod);
}

void NodeWiseMeshPartitioner::writeBinary(const std::string& file_name_base)
{
    writePropertiesBinary(file_name_base);
    const auto elem_integers = writeConfigDataBinary(file_name_base);

    const std::vector<IntegerType>& num_elem_integers
                                         = std::get<0>(elem_integers);
    const std::vector<IntegerType>& num_g_elem_integers
                                         = std::get<1>(elem_integers);
    writeElementsBinary(file_name_base, num_elem_integers,
                        num_g_elem_integers);

    writeNodesBinary(file_name_base);
}

void NodeWiseMeshPartitioner::writeConfigDataASCII
                                    (const std::string& file_name_base)
{
    const std::string fname =
                              file_name_base + "_partitioned_cfg"
                              + std::to_string(_npartitions) + ".msh";
    std::fstream os_subd_head(fname, std::ios::out | std::ios::trunc);
    const std::string mesh_info =
        "Subdomain mesh ("
        "Number of nodes; Number of base nodes;"
        " Number of regular elements; Number of ghost elements;"
        " Number of non-ghost base nodes; Number of non-ghost nodes"
        " Number of base nodes of the global mesh;"
        " Number of nodes of the global mesh;"
        " Number of integer variables to define non-ghost elements;"
        " Number of integer variables to define ghost elements.)";
    os_subd_head << mesh_info << "\n";
    os_subd_head << _npartitions << "\n";

    for (const auto& partition : _partitions)
    {
        os_subd_head << partition.nodes.size();
        os_subd_head << " " << partition.number_of_base_nodes;
        os_subd_head << " " << partition.regular_elements.size();
        os_subd_head << " " << partition.ghost_elements.size();
        os_subd_head << " " << partition.number_of_non_ghost_base_nodes;
        os_subd_head << " " << partition.number_of_non_ghost_nodes;
        os_subd_head << " " << _mesh->getNumberOfBaseNodes();
        os_subd_head << " " << _mesh->getNumberOfNodes();
        os_subd_head << " " << getNumberOfIntegerVariablesOfElements(
                                   partition.regular_elements);
        os_subd_head << " " << getNumberOfIntegerVariablesOfElements(
                                   partition.ghost_elements)
                     << " 0\n";
    }
}

void NodeWiseMeshPartitioner::writeElementsASCII(const std::string& file_name_base)
{
    const std::string fname = file_name_base + "_partitioned_elems_"
                              + std::to_string(_npartitions) + ".msh";
    std::fstream os_subd(fname, std::ios::out | std::ios::trunc);
    for (const auto& partition : _partitions)
    {
        // Set the local node indices of the current partition.
        IntegerType node_local_id_offset = 0;
        std::vector<IntegerType> nodes_local_ids(_mesh->getNumberOfNodes(), -1);
        for (const auto* node : partition.nodes)
        {
            nodes_local_ids[node->getID()] = node_local_id_offset;
            node_local_id_offset++;
        }

        for (const auto* elem : partition.regular_elements)
        {
            writeLocalElementNodeIndices(os_subd, *elem, nodes_local_ids);
        }
        for (const auto* elem : partition.ghost_elements)
        {
            writeLocalElementNodeIndices(os_subd, *elem, nodes_local_ids);
        }
        os_subd << std::endl;
    }
}

void NodeWiseMeshPartitioner::writeNodesASCII(const std::string& file_name_base)
{
    const std::string fname = file_name_base + "_partitioned_nodes_"
                              + std::to_string(_npartitions) + ".msh";
    std::fstream os_subd_node(fname, std::ios::out | std::ios::trunc);
    os_subd_node.precision(std::numeric_limits<double>::digits10);
    os_subd_node.setf(std::ios::scientific);

    for (const auto& partition : _partitions)
    {
        for (const auto* node : partition.nodes)
        {
            double const* coords = node->getCoords();
            os_subd_node << _nodes_global_ids[node->getID()] << " " << coords[0]
                         << " " << coords[1] << " " << coords[2] << "\n";
        }
        os_subd_node << std::endl;
    }
}

void NodeWiseMeshPartitioner::writeASCII(const std::string& file_name_base)
{
    writeConfigDataASCII(file_name_base);
    writeElementsASCII(file_name_base);
    writeNodesASCII(file_name_base);
}

void NodeWiseMeshPartitioner::getElementIntegerVariables(
    const MeshLib::Element& elem,
    const std::vector<IntegerType>& local_node_ids,
    std::vector<IntegerType>& elem_info,
    IntegerType& counter)
{
    unsigned mat_id = 0;  // TODO: Material ID to be set from the mesh data
    const IntegerType nn = elem.getNumberOfNodes();
    elem_info[counter++] = mat_id;
    elem_info[counter++] = static_cast<unsigned>(elem.getCellType());
    elem_info[counter++] = nn;

    for (IntegerType i = 0; i < nn; i++)
    {
        elem_info[counter++] = local_node_ids[elem.getNodeIndex(i)];
    }
}

void NodeWiseMeshPartitioner::writeLocalElementNodeIndices(
    std::ostream& os,
    const MeshLib::Element& elem,
    const std::vector<IntegerType>& local_node_ids)
{
    unsigned mat_id = 0;  // TODO: Material ID to be set from the mesh data
    os << mat_id << " " << static_cast<unsigned>(elem.getCellType())
       << " " << elem.getNumberOfNodes() << " ";
    for (unsigned i = 0; i < elem.getNumberOfNodes(); i++)
    {
        os << " " << local_node_ids[elem.getNodeIndex(i)];
    }
    os << "\n";
}

void NodeWiseMeshPartitioner::resetGlobalNodeIndices()
{
    for (std::size_t i = 0; i < _mesh->_nodes.size(); i++)
    {
        _mesh->_nodes[i]->setID(_nodes_global_ids[i]);
    }
    // sort
    std::sort(_mesh->_nodes.begin(), _mesh->_nodes.end(),
              [](const MeshLib::Node* a, const MeshLib::Node* b) {
                  return a->getID() < b->getID();
              });
}

void NodeWiseMeshPartitioner::writeGlobalMeshVTU(
    const std::string& file_name_base)
{
    resetGlobalNodeIndices();
    MeshLib::IO::VtuInterface writer(_mesh.get());
    const std::string npartitions_str = std::to_string(_npartitions);
    writer.writeToFile(file_name_base + "_node_id_renumbered_partitions_" +
                       npartitions_str + ".vtu");
}

}  // namespace MeshLib
