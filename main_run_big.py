
###########################
if __name__ == '__main__':
###########################
    # <rawcell>
    # define a specific fire to look at 
    # <codecell>


    reload(mt)
    reload(compute_transmittance)
    reload(call_models)
    reload(find_PBL)
    reload(msg)

    tagmakekmz=0
    tagplotraw=0 
    tagprint1=0  # Create raw figures of MIR-TIR images
    AOT_OMI_tag=0
    tagfigures=1 # Create figure os images

    #Study_case=sys.argv[1]
    #fire_name = str(Study_case)    # Input Parameters
    ###################
    myEmail = "ronan.paugam@kcl.ac.uk"
    fire_name  = 'SA11_07' # test


    flag_ecmwf = False

    # Southern Africa - August 2011 fire events
    if fire_name == 'SA11_01':
        create_out = True
        path_root      = '../data_SA11_01/'
        chernobyl_fire_lat =  -20.37
        chernobyl_fire_lon =  20.09
    if fire_name == 'SA11_02':
        create_out = True
        path_root      = '../data_SA11_02/'
        chernobyl_fire_lat =  -18.27
        chernobyl_fire_lon =  19.88
    if fire_name == 'SA11_03':
        create_out = True
        path_root      = '../data_SA11_03/'
        chernobyl_fire_lat =  -22.06
        chernobyl_fire_lon =  22.51
    if fire_name == 'SA11_04':
        create_out = True
        path_root      = '../data_SA11_04/'
        chernobyl_fire_lat =  -8.91
        chernobyl_fire_lon =  21.14
    if fire_name == 'SA11_05':
        create_out = True
        path_root      = '../data_SA11_05/'
        chernobyl_fire_lat =  -5.47
        chernobyl_fire_lon =  26.36

    if fire_name == 'SA11_06':
        create_out = True
        path_root      = '../data_SA11_06/'
        chernobyl_fire_lat =  -23.9
        chernobyl_fire_lon =  26.21
    if fire_name == 'SA11_07':
        create_out = True
        path_root      = '../data_SA11_07/'
        chernobyl_fire_lat =  -18.84
        chernobyl_fire_lon =  23.7
    if fire_name == 'SA11_08':
        create_out = True
        path_root      = '../data_SA11_08/'
        chernobyl_fire_lat =  -22.69
        chernobyl_fire_lon =  22.01
    if fire_name == 'SA11_09':
        create_out = True
        path_root      = '../data_SA11_09/'
        chernobyl_fire_lat =  -22.04
        chernobyl_fire_lon =  21.88
    if fire_name == 'SA11_10':
        create_out = True
        path_root      = '../data_SA11_10/'
        chernobyl_fire_lat =  -18.76
        chernobyl_fire_lon =  19.38
    if fire_name == 'SA11_11':
        create_out = True
        path_root      = '../data_SA11_11/'
        chernobyl_fire_lat =  -23.22
        chernobyl_fire_lon =  31.64
    if fire_name == 'SA11_12':
        create_out = True
        path_root      = '../data_SA11_12/'
        chernobyl_fire_lat =  -23.5
        chernobyl_fire_lon =  32.27
    if fire_name == 'SA11_13':
        create_out = True
        path_root      = '../data_SA11_13/'
        chernobyl_fire_lat =  -24.76
        chernobyl_fire_lon =  22.65
    if fire_name == 'SA11_14':
        create_out = True
        path_root      = '../data_SA11_14/'
        chernobyl_fire_lat =  -25.3
        chernobyl_fire_lon =  21.61
    if fire_name == 'SA11_15':
        create_out = True
        path_root      = '../data_SA11_15/'
        chernobyl_fire_lat =  -20.27
        chernobyl_fire_lon =  21.98
    if fire_name == 'SA11_16':
        create_out = True
        path_root      = '../data_SA11_16/'
        chernobyl_fire_lat =  -12.95
        chernobyl_fire_lon =  36.58
    if fire_name == 'SA11_17':
        create_out = True
        path_root      = '../data_SA11_17/'
        chernobyl_fire_lat =  -24.58
        chernobyl_fire_lon =  22.13
    if fire_name == 'SA11_18':
        create_out = True
        path_root      = '../data_SA11_18/'
        chernobyl_fire_lat =  -24.25
        chernobyl_fire_lon =  21.69
    if fire_name == 'SA11_19':
        create_out = True
        path_root      = '../data_SA11_19/'
        chernobyl_fire_lat =  -20.58
        chernobyl_fire_lon =  23.29
    if fire_name == 'SA11_20':
        create_out = True
        path_root      = '../data_SA11_20/'
        chernobyl_fire_lat =  -6.79
        chernobyl_fire_lon =  34.0
    if fire_name == 'SA11_21':
        create_out = True
        path_root      = '../data_SA11_21/'
        chernobyl_fire_lat =  -6.76
        chernobyl_fire_lon =  30.98
    if fire_name == 'SA11_22':
        create_out = True
        path_root      = '../data_SA11_22/'
        chernobyl_fire_lat =  -23.65
        chernobyl_fire_lon =  22.44
    if fire_name == 'SA11_23':
        create_out = True
        path_root      = '../data_SA11_23/'
        chernobyl_fire_lat =  -22.57
        chernobyl_fire_lon =  20.83
    if fire_name == 'SA11_24':
        create_out = True
        path_root      = '../data_SA11_24/'
        chernobyl_fire_lat =  -19.43
        chernobyl_fire_lon =  13.47


    if fire_name == 'SA11_25':
        create_out = True
        path_root      = '../data_SA11_25/'
        chernobyl_fire_lat =  -7.68
        chernobyl_fire_lon =  22.26
    if fire_name == 'SA11_26':
        create_out = True
        path_root      = '../data_SA11_26/'
        chernobyl_fire_lat =  -5.1
        chernobyl_fire_lon =  27.87
    if fire_name == 'SA11_27':
        create_out = True
        path_root      = '../data_SA11_27/'
        chernobyl_fire_lat =  -6.19
        chernobyl_fire_lon =  27.87
    if fire_name == 'SA11_28':
        create_out = True
        path_root      = '../data_SA11_28/'
        chernobyl_fire_lat =  -4.7
        chernobyl_fire_lon =  30.97
    if fire_name == 'SA11_29':
        create_out = True
        path_root      = '../data_SA11_29/'
        chernobyl_fire_lat =  -18.77
        chernobyl_fire_lon =  25.65
    if fire_name == 'SA11_30':
        create_out = True
        path_root      = '../data_SA11_30/'
        chernobyl_fire_lat =  -18.77
        chernobyl_fire_lon =  25.65
    if fire_name == 'SA11_31':
        create_out = True
        path_root      = '../data_SA11_31/'
        chernobyl_fire_lat =  -16.86
        chernobyl_fire_lon =  18.73


    # Ronans example
    if fire_name == 'TT15':
        create_out = True
        path_root      = '../data_test/'
        chernobyl_fire_lat =  51.27
        chernobyl_fire_lon =  29.8173    


    #test ecacces certificate
    if flag_ecmwf:
        ecmwf_tools.test_certificate()


    #define area for plot
    #N_pixel = 400 # number of pixels seletected in the granule around the fire location
    N_pixel = 400 # number of pixels seletected in the granule around the fire location
    #zoom_size = 40  # number of pixel in the zoom box ploting the raw data
    zoom_size = 100  # number of pixel in the zoom box ploting the raw data


    #define area to georeference
    dxy                 = 500. # resolution of the georeferenced ROI
    #georef_roi_box_size = 60   # extension of the ROI box on each side of the fire
    georef_roi_box_size = 120   # extension of the ROI box on each side of the fire

    windowsize=georef_roi_box_size*2*(1000/dxy)

    #path directory. Here they are all link to current dir.
    path_MOYDIS03  = path_root+'MODIS/'
    path_MOYDIS021 = path_root+'MODIS/'
    path_MOYDIS14  = path_root+'MODIS/'
    path_MOYDIS04  = path_root+'MODIS/'
    path_ECMWF  = path_root+'ECMWF/'
    path_OMAERO  = path_root+'OMAERO/'
   
    if flag_ecmwf:
        ensure_dir(path_ECMWF)
        ensure_dir(path_ECMWF+'Analysis/')
        ensure_dir(path_ECMWF+'Atmospheric_profile/')
        ensure_dir(path_ECMWF+'HbL/')

    path_out       = path_root+'Output/'+fire_name+'/'
    path_image     = path_root+'Output/'+fire_name+'/images/'
    path_kmz       = path_root+'Output/'+fire_name+'/KMZ/'
    

    #create directory if not present
    mt.ensure_dir(path_MOYDIS03)
    mt.ensure_dir(path_MOYDIS021)
    mt.ensure_dir(path_MOYDIS14)
    mt.ensure_dir(path_MOYDIS04)
    mt.ensure_dir(path_OMAERO)
    mt.ensure_dir(path_out)
    mt.ensure_dir(path_image)
    mt.ensure_dir(path_kmz)

    sat= []; year= []; doy= []; dateTime= []; time= []; fire_lat= []; fire_lon = []
    
    # overpass info is coming from http://cloudsgate2.larc.nasa.gov/cgi-bin/predict/predict.cgi
    data_terra = load_overpass_info(path_root + 'Input/'+fire_name+'_Terra_overpasses.txt')
    data_terra_frame = pandas.DataFrame(data_terra)
    
    data_aqua  = load_overpass_info(path_root + 'Input/'+fire_name+'_Aqua_overpasses.txt')
    data_aqua_frame = pandas.DataFrame(data_aqua)
    data_aqua_frame.index = np.arange(data_terra_frame.index[-1]+1,data_terra_frame.index[-1]+1+len(data_aqua_frame.index))
    data_all = pandas.concat([data_terra_frame, data_aqua_frame])

    sat_name = len(data_terra)*['T'] + len(data_aqua)*['A']
    for i in range(data_all.shape[0]):
	print i
        sat.append(sat_name[i])
        year.append(int(data_all['col1'][i]))
        month_  = int(data_all['col2'][i])
        day_    = int(data_all['col3'][i])
        hour_   = int(data_all['col4'][i])
        minute_ = int(data_all['col5'][i])
	if data_all['col6'][i]==60:
             data_all['col6'][i]=59
        second_ = int(data_all['col6'][i])
        date_ = datetime.datetime(year[-1],month_,day_,hour_,minute_,second_)
        dateTime.append(date_)
        #doy.append(int(date_.strftime('%j')))
        #time.append(int(str(data_all['col4'][i])+str(data_all['col5'][i])))
        fire_lat.append(chernobyl_fire_lat)
        fire_lon.append(chernobyl_fire_lon)


    #upload MSG FRP data
    #######################################

    fire_tab=msg.upload_msg(data_all, georef_roi_box_size, chernobyl_fire_lat, chernobyl_fire_lon)

    #MSG_Fire_array=msg.read_msg(Fire_MSG,year,month,day,time,conv_ll2utm)


    # loop over all overpasses
    ##########################

    if create_out :
        if flag_ecmwf:
            #fetch PBL data from ECMWF server.
            obs_hbl = find_PBL.upload_pbl_from_ecmwf(fire_name,sat,dateTime,fire_lat,fire_lon,path_ECMWF,path_out)
        else:
            obs_hbl = None

	AOT_images=np.zeros((windowsize, windowsize, data_all.shape[0]))
	ImgCov=np.zeros((data_all.shape[0], 3))
        item = tuple(['XX']+[datetime.datetime.now()]+['XX','XX','XX'] + [0.,]*36 + ['XX','XX'] + [0.]*6)
        out = np.array([item]*data_all.shape[0],dtype=np.dtype([
                                                ('fireName','S100'),('fireDatetime',datetime.datetime),('fireDate','S100'),('fireTime','S100'),('satellite','S1'),  \
                                                ('fireLon_inF',float),('fireLat_inF',float),('viewAngle',float),('azimuthAngle',float),                             \
                                                ('FRP_kcl_iF_MW',float), ('FRP_mod_iF_MW',float),('afArea_iF_ha',float),('BT_iF_K',float),                          \
                                                ('FRP_kcl_bF_MW',float), ('FRP_mod_bF_MW',float),('afArea_bF_ha',float),('BT_bF_K',float),                          \
                                                ('FRP_kcl_inF_MW',float),('FRP_mod_inF_MW',float),('afArea_inF_ha',float),('BT_inF_K',float),                       \
                                                ('FRP_kcl_bnF_MW',float), ('FRP_mod_bnF_MW',float),('afArea_bnF_ha',float),('BT_bnF_K',float),                      \
                                                ('iF_L4',float),('iF_L4b',float),('iF_L11',float),('iF_L11b',float),                                                \
                                                ('bF_L4',float),('bF_L4b',float),('bF_L11',float),('bF_L11b',float),                                                \
                                                ('inF_L4',float),('inF_L4b',float),('inF_L11',float),('inF_L11b',float),                                            \
                                                ('bnF_L4',float),('bnF_L4b',float),('bnF_L11',float),('bnF_L11b',float),                                            \
                                                ('path_geoRef_MIR','S200'),('path_geoRef_ColorComp','S200'),                                                        \
                                                ('InJH_prmv0',float),('InJH_prmv1',float),('InJH_prmv2',float),('InJH_sof',float),                                  \
                                                ('hbl_ecmwf',float),('hbl_kcl',float)                                                                             ]))
        out = out.view(np.recarray)


        items = tuple([0.0]*18)
        MODIS_AOT_array = np.array([items]*data_all.shape[0],dtype=np.dtype([('N_AOT',float),('Sum_AOT',float),('P1_AOT',float),\
		('P5_AOT',float),('P10_AOT',float),('P50_AOT',float),('P90_AOT',float),('P95_AOT',float),('P99_AOT',float),\
		('N_AOTuN',float),('Sum_AOTuN',float),('P1_AOTuN',float),('P5_AOTuN',float),('P10_AOTuN',float),('P50_AOTuN',float),\
		('P90_AOTuN',float),('P95_AOTuN',float),('P99_AOTuN',float)]))   


        MODIS_AOT_array = MODIS_AOT_array.view(np.recarray)


        for ii, [sat_, dateTime_, fire_lat_, fire_lon_] in enumerate(zip(sat,dateTime,fire_lat,fire_lon)):
            print ''
            print ii, dateTime_.year, dateTime_.month, dateTime_.day, '@' , '{:02d}:{:02d}'.format(dateTime_.hour,dateTime_.minute)
            print '###############'
            
            out_ , grid_e, grid_n, dxy, conv_ll2utm, conv_utm2ll, AIcha, AIUncha, AISta, AIUnSta, Aero_geo, percent_obs= load_georef_plot(fire_name, sat_,\
							  dateTime_, fire_lat_,fire_lon_,obs_hbl,\
                                                          flag_ecmwf,path_ECMWF,path_image, fire_tab)

            AOT_images[:,:,ii]=Aero_geo[:,:]
            ImgCov[ii,:]=percent_obs


            MODIS_AOT_array.N_AOT[ii]=AIcha[0]
            MODIS_AOT_array.Sum_AOT[ii]=AIcha[1]
            MODIS_AOT_array.P1_AOT[ii]=AISta[0]
            MODIS_AOT_array.P5_AOT[ii]=AISta[1]
            MODIS_AOT_array.P10_AOT[ii]=AISta[2]
            MODIS_AOT_array.P50_AOT[ii]=AISta[3]
            MODIS_AOT_array.P90_AOT[ii]=AISta[4]
            MODIS_AOT_array.P95_AOT[ii]=AISta[5]
            MODIS_AOT_array.P99_AOT[ii]=AISta[6]
            MODIS_AOT_array.N_AOTuN[ii]=AIUncha[0]
            MODIS_AOT_array.Sum_AOTuN[ii]=AIUncha[1]
            MODIS_AOT_array.P1_AOTuN[ii]=AIUnSta[0]
            MODIS_AOT_array.P5_AOTuN[ii]=AIUnSta[1]
            MODIS_AOT_array.P10_AOTuN[ii]=AIUnSta[2]
            MODIS_AOT_array.P50_AOTuN[ii]=AIUnSta[3]
            MODIS_AOT_array.P90_AOTuN[ii]=AIUnSta[4]
            MODIS_AOT_array.P95_AOTuN[ii]=AIUnSta[5]
            MODIS_AOT_array.P99_AOTuN[ii]=AIUnSta[6]

# need to put something for processing the AOT
            
            for field in out.dtype.names:
                if field == 'satellite':
                    out[field][ii]= sat_ 
                else:
                    out[field][ii] = out_[field][0]

	    #if ii==data_all.shape[0]-1
		
        #MSG_Fire_array=msg.read_msg(Fire_MSG,year,month,day,time,conv_ll2utm)
        MSG_Fire_array=msg.read_msg_ts(fire_tab,conv_ll2utm, grid_e, grid_n)
        FRE, FREun, conta = msg.plot_dados(MSG_Fire_array, sat, dateTime, path_image, fire_name, 1)
        LCdom_geo, LCdom_per_geo = msg.plot_fire_lc(MSG_Fire_array, dateTime, fire_lat[0], fire_lon[0], N_pixel, conv_ll2utm, conv_utm2ll, grid_n, grid_e, dxy, path_image, fire_name)
        percent_fire_count, percent_fRP_count, classes, dias, lcp = msg.plot_hist_frplc(MSG_Fire_array, grid_e, grid_n, LCdom_geo, path_image, fire_name)
        msg.plot_aot_frp(MODIS_AOT_array, ImgCov, FRE, FREun, dateTime, sat, path_image, fire_name)
        msg.pltsct_aot_frp(MODIS_AOT_array, ImgCov, FRE, FREun, dateTime, sat, path_image, fire_name)
        msg.plot_aot_prop(MODIS_AOT_array, ImgCov, dateTime, sat, path_image, fire_name)
        #remove empty obs. (does not mean no fire, but nothing at all.)
        idx = np.where(out.path_geoRef_ColorComp != 'XX') # obs where there was nothing. outisde modis track
        out = out[idx]
        out = out.view(np.recarray)

        #sort obs per time
        out.sort(order='fireDatetime')
        out_4dump = drop_fields(out, 'fireDatetime')
        out_4dump = out_4dump.view(np.recarray)

        MSG_Fire_array.sort(order='date')
        MSG_Fire_array_simple=drop_fields(MSG_Fire_array, 'date')
        MSG_Fire_array_simple=MSG_Fire_array_simple.view(np.recarray)

        MODIS_AOT_array=MODIS_AOT_array.view(np.recarray)
        classes=classes.view(np.recarray)
        #dump fire info text format
        df = pandas.DataFrame.from_records(out_4dump)
        df.to_csv(path_out+fire_name+'fire_info_all.csv',sep=',')

        #dump npy format
        np.save(path_out+fire_name+'fire_info_all',[out, grid_e, grid_n, dxy])
        np.save(path_out+fire_name+'fire_info_MSG',[MSG_Fire_array, FRE, FREun, conta, lcp])
        np.save(path_out+fire_name+'AOT_info_MSG',[AOT_images, MODIS_AOT_array])
        np.save(path_out+fire_name+'LC_info_MSG',[percent_fire_count, percent_fRP_count, dias, classes, LCdom_geo, LCdom_per_geo])


        #save data into matlab files
        scipy.io.savemat(path_out+'Data_geo_'+fire_name+'.mat', {'grid_e':grid_e, 'grid_n':grid_n, 'LCdom_geo':LCdom_geo, 'LCdom_per_geo':LCdom_per_geo})
        scipy.io.savemat(path_out+'Data_Fire_'+fire_name+'.mat', {'MSG_Fire_array_simple':MSG_Fire_array_simple, 'Out':out_4dump, 'FRE':FRE, 'FREun':FREun, 'conta':conta, 'LCf':lcp})
        scipy.io.savemat(path_out+'Data_AOT_'+fire_name+'.mat', {'AOT_images':AOT_images, 'MODIS_AOT_array':MODIS_AOT_array})
        scipy.io.savemat(path_out+'LCov_'+fire_name+'.mat', {'Per_Counts':percent_fire_count, 'Per_FRP':percent_fRP_count, 'classes':classes, 'dias':dias})
	

        #percent_fire_count, percent_fRP_count, classes, dias, lcp

    
    else:
        out, grid_e, grid_n, dxy = np.load(path_out+fire_name+'fire_info_all.npy')
        out = out.view(np.recarray)

