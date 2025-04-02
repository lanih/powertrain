function voltage = battery_model(current, p_count, s_count) % eventually add temperature, SoC
    
    cell_voltage = [4.2028
    4.1471
    4.1020
    4.0297
    3.9603
    3.8525
    3.7800
    3.7337
    3.7003
    3.6667
    3.6201
    3.5473
    3.4748
    3.3018
    3.0093
    2.9342]


    surfp = [0.4240
    0.4458
    0.4652
    0.5041
    0.5429
    0.6007
    0.6584
    0.7162
    0.7739
    0.8317
    0.8894
    0.9472
    0.9706
    0.9847
    0.9988
    1.0017]

    surfn = [0.9867
    1.0023
    1.0159
    1.0429
    1.0699
    1.1101
    1.1503
    1.1905
    1.2307
    1.2709
    1.3111
    1.3513
    1.3676
    1.3774
    1.3872
    1.3892]
    
    pack_voltage = cell_voltage * s_count; % Eventually use to model battery drain 

    

