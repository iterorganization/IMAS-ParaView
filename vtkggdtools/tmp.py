from benchmarks.benchmarks import Convert, ConvertHDF5, ConvertMDSPlus

converter = Convert()
converter.setup()
converter.time_convert()

converthdf5 = ConvertHDF5()
converthdf5.setup()
converthdf5.time_load_and_convert()
converthdf5.time_load_and_convert_lazy()

convertmdsplus = ConvertMDSPlus()
convertmdsplus.setup()
convertmdsplus.time_load_and_convert()
convertmdsplus.time_load_and_convert_lazy()
