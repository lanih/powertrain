function pack_voltage = BatteryPack (cell, series_cells, parallel_cells)

series_cells = 16; 
parallel_cells = 3; 

pack_voltage = cell.volt(9)*series_cells


end