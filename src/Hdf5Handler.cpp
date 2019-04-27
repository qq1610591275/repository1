#ifdef FORMAT_HDF5

#include "Hdf5Handler.h"

#include "H5Cpp.h"
#include "hdf5.h"
#include "hdf5_hl.h"

#include "log.h"
#include "Material.h"

namespace antmoc
{

Hdf5Handler::Hdf5Handler() {}

Hdf5Handler::~Hdf5Handler() {}


/**
 * @brief 获取材料数组
 * @retrun 返回材料数组
 */
std::map<int, Material *> Hdf5Handler::getMaterials() {
	return _materials;
}


/**
 * @breif: 从H5文件读取宏观材料数据
 */
std::map<int, Material *> Hdf5Handler::readMGXSFromFile(std::string filename)
{int count =0 ;
  H5::H5File *h5_file = nullptr;
  std::stringstream intstr;
  hsize_t dims[2];
  int row = 0, j = 0, energy_groups = 0;
  std::string groupeName = "", datasetName = "";
  // Open the file
  h5_file = new H5::H5File(filename, H5F_ACC_RDONLY);
  // Get the number of energy groups
  H5::Attribute attr_energy = h5_file->openAttribute("energy_groups");
  H5::IntType attr_energy_type = attr_energy.getIntType();
  attr_energy.read(attr_energy_type, &energy_groups);

  hid_t num_groups = h5_file->getNumObjs();
  for (int i = 1; i <= num_groups; i++)
  {
    intstr.str("");
    if (i < 10)
    {
      intstr << i;
      groupeName = "  " + intstr.str();
    }
    else if (i < 100)
    {
      intstr << i;
      groupeName = " " + intstr.str();
    }
    else
    {
      intstr << i;
      groupeName = intstr.str();
    }
    if (!h5_file->nameExists(groupeName))
      log_printf(ERROR, "the currrent group name %s doesn't exsit!", groupeName);
    H5::Group group = h5_file->openGroup(groupeName);
    hid_t num_datasets = group.getNumObjs();
    log_printf(DEBUG, "the half dataset numbers is %d", num_datasets / 2);
    if (group.nameExists("  0-reactions"))
      num_datasets = num_datasets / 2 - 1;
    else
      num_datasets = num_datasets / 2;
	int temp = num_datasets,count_num=0; 
	while(temp>10){   // calculate the number of digits
		count_num++;
		temp /= 10;
	}
	count_num++;
    for (j = 1; j <= num_datasets; j++)
    {
      intstr.str("");
      std::string dataName = "", scatterName = "";
      if (j < 10)
      {
        intstr << j;
        dataName = "  " + intstr.str();
      }
      else if (j < 100)
      {
        intstr << j;
        dataName = " " + intstr.str();
      }
      else
      {
        intstr << j;
        dataName = intstr.str();
      }
      scatterName = dataName + "-scattering";
      dataName = dataName + "-reactions";
      int mater_id = 0;
      mater_id = i*pow(10,count_num)+j; // transform the string to integer:1-1->11

      // read the reactions data from the current group
      H5::DataSet dataset_react(group.openDataSet(dataName));
      H5::DataSpace filespace_react = dataset_react.getSpace();
      hid_t rank = filespace_react.getSimpleExtentDims(dims);
      H5::DataSpace mspace_react(rank, dims);
      H5::FloatType floattype_react(dataset_react);
      double *value_data = (double *)malloc(dims[0] * dims[1] * sizeof(double));
      if (dims[0] != energy_groups)
        log_printf(ERROR, "the format of data is incorrect");
      double sig_a[dims[0]] = {0.0}, sig_c[dims[0]] = {0.0}, sig_f[dims[0]] = {0.0}, sig_t[dims[0]] = {0.0}, sig_nu_f[dims[0]] = {0.0};
      dataset_react.read(value_data, floattype_react, mspace_react, filespace_react);
      for (row = 0; row < dims[0]; row++)
      {
        sig_a[row] = value_data[row * dims[1] + 0];
        sig_c[row] = value_data[row * dims[1] + 1];
        sig_f[row] = value_data[row * dims[1] + 2];
        sig_nu_f[row] = value_data[row * dims[1] + 3];
        sig_t[row] = value_data[row * dims[1] + 0];
      }
      delete [] value_data; // free the mem

      // read the scatter date from current group
      H5::DataSet dataset_sca(group.openDataSet(scatterName));
      H5::DataSpace filespace = dataset_sca.getSpace();
      rank = filespace.getSimpleExtentDims(dims);
      H5::DataSpace mspace1(rank, dims);
      H5::FloatType floattype(dataset_sca);
      double array[dims[0] * dims[1]] = {0.0};
      if (dims[0] != energy_groups || dims[1] != energy_groups)
        log_printf(ERROR, "the data format is incorrect");
      dataset_sca.read(array, floattype, mspace1, filespace);
      for (row = 0; row < dims[0]; ++row)
      {
        if ((dims[1] * row + 5) % 5 == 0)
          sig_t[row] += array[dims[1] * row + 5];
      }
      // store the info to the new material
	  count++;
      Material *material = new Material(mater_id);
      material->setNumEnergyGroups(energy_groups);
      material->setSigmaT(sig_t, energy_groups);
      material->setNuSigmaF(sig_nu_f, energy_groups);
      material->setSigmaF(sig_f, energy_groups);
      material->setChi(sig_c, energy_groups);
      material->setSigmaS(array, energy_groups *energy_groups);
      _materials.insert({material->getId(), material});
    }
  }
  return _materials;
}


/**
  * @brief 通过材料的id获取材料对象
  * @parm  mat_id是材料的id
  * @return 返回找到的材料
  */
Material *Hdf5Handler::findMaterial(int mat_id) {
	// 第二个人添加的注释
	if (_materials.find(mat_id) == _materials.end())
		log_printf(ERROR, "can't find the material that id is %d", mat_id);
	return _materials.at(mat_id);
}


} /* namespace antmoc */

#endif /* FORMAT_HDF5 */
