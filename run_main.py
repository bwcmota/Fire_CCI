
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
    #fire_name  = 'SA11_05' # test


    flag_ecmwf = False

    # Southern Africa - August 2011 fire events
    fire_name = 'SA11_15'
    create_out = True
    path_root      = '../data_SA11_15/'
    chernobyl_fire_lat =  -25.3-0.0 #1.75
    chernobyl_fire_lon =  21.61+0.5
    sat='A'
    ano_i=2011
    mes_i=10
    dia_i=15
    hora_i=12
    ie=48

    tagprocess=1

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
    georef_roi_box_size = 160   # extension of the ROI box on each side of the fire

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

    sat= []; year= []; doy= []; dateTime= []; time= []; fire_lat= []; fire_lon = []; datedelta=[]; vez=[]
    
    # overpass info is coming from http://cloudsgate2.larc.nasa.gov/cgi-bin/predict/predict.cgi
    data_terra = load_overpass_info(path_root + 'Input/'+fire_name+'_Terra_overpasses.txt')
    data_terra_frame = pandas.DataFrame(data_terra)
    
    data_aqua  = load_overpass_info(path_root + 'Input/'+fire_name+'_Aqua_overpasses.txt')
    data_aqua_frame = pandas.DataFrame(data_aqua)
    data_aqua_frame.index = np.arange(data_terra_frame.index[-1]+1,data_terra_frame.index[-1]+1+len(data_aqua_frame.index))
    data_all = pandas.concat([data_terra_frame, data_aqua_frame])
    sat_name = len(data_terra)*['T'] + len(data_aqua)*['A']

    timetag=datetime.datetime(ano_i, mes_i, dia_i, hora_i)
    for i in range(data_all.shape[0]):
	print i
	vez.append(i)
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
        delta_=date_-timetag
        dateTime.append(date_)
        datedelta.append(delta_)
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


#        for ii, [sat_, dateTime_, fire_lat_, fire_lon_] in enumerate(zip(sat,dateTime,fire_lat,fire_lon)):
#            print ''
#            print ii, dateTime_.year, dateTime_.month, dateTime_.day, '@' , '{:02d}:{:02d}'.format(dateTime_.hour,dateTime_.minute)
#            print '###############'

        datedelta=abs(np.array(datedelta))
        vez=np.array(vez)
        flag=vez[datedelta==datedelta.min()][0]

        sat_=sat[flag]
        dateTime_=dateTime[flag]
        fire_lat_=chernobyl_fire_lat
        fire_lon_=chernobyl_fire_lon 

        out_ , grid_e, grid_n, dxy, conv_ll2utm, conv_utm2ll, AIcha, AIUncha, AISta, AIUnSta, Aero_geo, percent_obs= load_georef_plot(fire_name, sat_,\
							  dateTime_, fire_lat_,fire_lon_,obs_hbl,\
                                                          flag_ecmwf,path_ECMWF,path_image, fire_tab)

	if tagprocess==1:

		print "REPLOJECT MSG FIRES"

		MSG_Fire_array=msg.read_msg_ts(fire_tab,conv_ll2utm, grid_e, grid_n, dateTime)
	#        FRE, FREun, conta = msg.plot_dados(MSG_Fire_array, sat, dateTime, path_image, fire_name, 1)
		print "GET ENERGY"
		FRE, FREun, conta, resetime = msg.cal_FRE(MSG_Fire_array, sat, dateTime)
		print "GET LANDCOVER MAPS"

		LCdom_geo, LCdom_per_geo = msg.cal_fire_lc(MSG_Fire_array, dateTime, fire_lat[0], fire_lon[0], N_pixel, conv_ll2utm, conv_utm2ll, grid_n, grid_e, dxy)
	#        percent_fire_count, percent_fRP_count, classes, dias, lcp = msg.plot_hist_frplc(MSG_Fire_array, grid_e, grid_n, LCdom_geo, path_image, fire_name)
		print "GET LANDCOVER STATS"
		percent_fire_count, percent_fRP_count, classes, dias, lcp, countlc = msg.cal_hist_frplc(MSG_Fire_array, grid_e, grid_n, LCdom_geo, resetime)

		perc_case=countlc[flag,:]/sum(countlc[flag,:])*100.0
		fre_case=FRE[flag]
		freun_case=FREun[flag]
		aot_case=AIcha
		aots_case=AISta
		aotun_case=AIUncha
		aotuns_case=AIUnSta

		scipy.io.savemat(path_out+'Case_Data_geo_'+str(ie)+'.mat', {'perc_case':perc_case, 'fre_case':fre_case, 'freun_case':freun_case, 'aot_case':aot_case, 'aots_case':aots_case, 'aotun_case':aotun_case, 'aotuns_case':aotuns_case})


		print str(fre_case), str((aot_case[1]-(aot_case[0]*aots_case[2]))[0])



