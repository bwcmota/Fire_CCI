
##########################################################
def load_georef_plot(fire_name,sat,date_fire,fire_lat,fire_lon,obs_hbl,flag_ecmwf,ecmwf_dir,path_image, fire_tab):


   #out_ , grid_e, grid_n, dxy = load_georef_plot(fire_name, sat_, dateTime_, fire_lat_,fire_lon_,obs_hbl,\
   #                                                       flag_ecmwf,path_ECMWF,path_image)

    #sat=sat_
    #date_fire=dateTime_


    # <rawcell>
    ############################################################################################
    # Finish intialisation of variable and direactory paths
    ############################################################################################
    # <codecell>

    #compute day and month from doy
    year = date_fire.year
    doy  = int(date_fire.strftime('%j'))
    time = 100*date_fire.hour + date_fire.minute   
    #modis time
    date_modis = datetime.datetime.strptime('{:04d}  {:03d} {:02d} {:02d}'.format(year,doy,int(time/100),5*(int((time-100*int(time/100)))/5)), '%Y %j %H %M') 
    time = 100*date_modis.hour + date_modis.minute
    month= date_modis.month
    day  = date_modis.day

    percent_obs = np.zeros(3)


    if sat[0] != 'A' : 
        sat_     = 'Terra'
        prefix21 = 'MOD021KM'
        prefix03 = 'MOD03'
        prefix14 = 'MOD14'
        prefix04 = 'MOD04_L2'
        dir_nasa03 = 'MOD03'
        dir_nasa21 = 'MOD021KM'
        dir_nasa14 = 'MOD14'
        dir_nasa04 = 'MOD04_L2'
    else: 
        sat_     = 'Aqua'
        prefix21 = 'MYD021KM'
        prefix03 = 'MYD03'
        prefix14 = 'MYD14'
        prefix04 = 'MYD04_L2'
        prefixOmaero = 'OMAEROG.003'
        dir_nasa03 = 'MYD03'
        dir_nasa21 = 'MYD021KM'
        dir_nasa14 = 'MYD14'
        dir_nasa04 = 'MYD04_L2'
    	dir_OMAERO = 'OMAEROG.003'

    # Set projection type and conversion tools      
    # tool to convert lon lat to UTM
    wgs84 = osr.SpatialReference( ) # Define a SpatialReference object
    wgs84.ImportFromEPSG( 4326 ) # And set it to WGS84 using the EPSG code
    utm = osr.SpatialReference()
    utm.SetWellKnownGeogCS( 'WGS84' )
    utm.SetUTM(mt.UTMZone(fire_lon, fire_lat))
    conv_ll2utm = osr.CoordinateTransformation(wgs84, utm)
    conv_utm2ll = osr.CoordinateTransformation(utm,wgs84) 

    AIcha=np.empty(2)
    AIcha[:]=np.NAN
    AIUncha=np.empty(2)
    AIUncha[:]=np.NAN
    AISta=np.empty(7)
    AISta[:]=np.NAN
    AIUnSta=np.empty(7)
    AIUnSta[:]=np.NAN
    windowsize=georef_roi_box_size*2*(1000/dxy)
    Aero_geo=np.empty((windowsize,windowsize))
    Aero_geo[:,:]=np.NAN
    # <headingcell level=1>
    ############################################################################################
    # 1. Upload the MOD03, MOD021km, and MOD14 MODIS product
    ############################################################################################
    # <codecell>

    ###############################################
    #upload MOD03 - MODIS geolocation files
    ###############################################
    MOD03_tmp_T  = prefix03 + '.A{:04d}{:03d}.{:04d}.005*.hdf'.format(year,doy,time)
    res = glob.glob( path_MOYDIS03 + MOD03_tmp_T)
    if len(res) == 0 :
        ftp = ftplib.FTP("ladsweb.nascom.nasa.gov")
        ftp.login("anonymous", myEmail)
        res = mt.download_atNasa(ftp,MOD03_tmp_T,year, doy, dir_nasa03, path_MOYDIS03)
        ftp.quit()
    MOD03 = res[0]    

    #load lat lon from MOD03
    gdal_data_lat = gdal.Open ('HDF4_SDS:UNKNOWN:'+MOD03+':0')
    lat03_granule = gdal_data_lat.ReadAsArray()
    gdal_data_lon = gdal.Open ('HDF4_SDS:UNKNOWN:'+MOD03+':1')
    lon03_granule = gdal_data_lon.ReadAsArray()
    Nl03, Ns03 =  lon03_granule.shape # Nbre line and sample

    #and view angle and azimuth
    gdal_data_azimuth = gdal.Open ('HDF4_SDS:UNKNOWN:'+MOD03+':4')
    azimuth = gdal_data_azimuth.ReadAsArray()*float(gdal_data_azimuth.GetMetadata()['scale_factor'])
    gdal_data_zenith = gdal.Open ('HDF4_SDS:UNKNOWN:'+MOD03+':3')
    zenith = gdal_data_zenith.ReadAsArray()*float(gdal_data_zenith.GetMetadata()['scale_factor'])
    

    ###############################################
    #upload MOD021 - reflectance data
    ###############################################
    MOD021_tmp_T = prefix21 +  '.A{:04d}{:03d}.{:04d}.005*.hdf'.format(year,doy,time)
    res = glob.glob( path_MOYDIS021 + MOD021_tmp_T)
    if len(res) == 0 :
        ftp = ftplib.FTP("ladsweb.nascom.nasa.gov")
        ftp.login("anonymous", "ronan.paugam@kcl.ac.uk")
        res = mt.download_atNasa(ftp,MOD021_tmp_T,year, doy, dir_nasa21,path_MOYDIS021)
        ftp.quit()        
    MOD021 = res[0]    

    # load TIR Brightness Temperature (BT) (band31) 
    print 'load TIR band'
    gdal_data = gdal.Open ('HDF4_EOS:EOS_SWATH:'+MOD021+':MODIS_SWATH_Type_L1B:EV_1KM_Emissive')
    iband = 10
    L_tir = np.array(gdal_data.ReadAsArray()[iband,:,:],dtype=np.float) * float(gdal_data.GetMetadata()['radiance_scales'].split(',')[iband])
    #apply transmittance
    L_tir = L_tir 
    #convert to brightness temperature
    wavelength = np.zeros(L_tir.shape) + 11.009
    BT_tir = mt.planck_temperature(wavelength,L_tir)
    # and apply a inverse of a model spectral function directly on BT to get central BT
    BT_tir = mt.apply_inv_SRF_model_on_BT(BT_tir,'TIR')

    #and MIR BT (band21)
    print 'load MIR band'
    iband = 1
    L_mir = np.array(gdal_data.ReadAsArray()[iband,:,:],dtype=np.float) * float(gdal_data.GetMetadata()['radiance_scales'].split(',')[iband])
    #apply transmittance
    L_mir = L_mir 
    #convert to brightness temperature
    wavelength = np.zeros(L_mir.shape) + 3.959
    BT_mir = mt.planck_temperature(wavelength,L_mir)
    # and apply the inverse of a model spectral function directly on BT to get central BT
    BT_mir = mt.apply_inv_SRF_model_on_BT(BT_mir,'MIR')


    # for true color composite imageload band1(red, 0.6 micron), 3(blue, 0.47 micron), 4(green, 0.55 micron)
    print 'load true color composite band'
    gdal_data = gdal.Open ('HDF4_EOS:EOS_SWATH:'+MOD021+':MODIS_SWATH_Type_L1B:EV_250_Aggr1km_RefSB')
    iband = 0
    band1 = np.array(gdal_data.ReadAsArray()[iband,:,:],dtype=np.float) * float(gdal_data.GetMetadata()['reflectance_scales'].split(',')[iband])
    gdal_data = gdal.Open ('HDF4_EOS:EOS_SWATH:'+MOD021+':MODIS_SWATH_Type_L1B:EV_500_Aggr1km_RefSB')
    iband = 0
    band3 = np.array(gdal_data.ReadAsArray()[iband,:,:],dtype=np.float)* float(gdal_data.GetMetadata()['reflectance_scales'].split(',')[iband])
    iband = 1
    band4 = np.array(gdal_data.ReadAsArray()[iband,:,:],dtype=np.float)* float(gdal_data.GetMetadata()['reflectance_scales'].split(',')[iband])
    size_im = [band1.shape[0],band1.shape[1],3] 
    imColorComp = np.zeros(size_im)#, dtype=np.uint8)
    imColorComp[:,:,0] = band1
    imColorComp[:,:,1] = band4
    imColorComp[:,:,2] = band3

    '''
    #get view angle of the fire location
    dist_ = np.abs(lat03_granule-fire_lat)**2 + np.abs(lon03_granule-fire_lon)**2
    min_dist = dist_.min()
    idx_min = np.where(dist_==min_dist)
    mod_pix_size,vzen,scananglemap = mt.modisImageInfoModel()
    fire_viewAngle = vzen.T[idx_min][0]
    '''

    #####################################################
    #upload the MOD14 product and filter MSG FRP product
    #####################################################
    MOD14_tmp_T  = prefix14 + '.A{:04d}{:03d}.{:04d}.005*.hdf'.format(year,doy,time)
    res = glob.glob( path_MOYDIS14 + MOD14_tmp_T)
    if len(res) == 0 :
        ftp = ftplib.FTP("ladsweb.nascom.nasa.gov")
        ftp.login("anonymous", myEmail)
        res = mt.download_atNasa(ftp,MOD14_tmp_T,year, doy, dir_nasa14, path_MOYDIS14)
        ftp.quit()
    MOD14 = res[0] 

    #load data
    print 'load MOD14'
    treshold_fp_detection = 3  # default value is 3
    fire_product = mt.read_mod14(MOD14,                             \
                                 year,month,day,time,               \
                                 treshold_fp_detection, conv_ll2utm )

    #Filter the MSG FRP data up to satellite overpass
    fire_MSG = msg.read_msg(fire_tab, year,month,day,time, conv_ll2utm)





    ######################################## 
    #upload MOD04
    ########################################    
    MOD04_L2_tmp_T = prefix04 +  '.A{:04d}{:03d}.{:04d}.006*.hdf'.format(year,doy,time)
    res = glob.glob( path_MOYDIS04 + MOD04_L2_tmp_T)
    if len(res) == 0 :
        ftp = ftplib.FTP("ladsweb.nascom.nasa.gov")
        ftp.login("anonymous", "ronan.paugam@kcl.ac.uk")
        res = mt.download_atNasa6(ftp,MOD04_L2_tmp_T,year, doy, dir_nasa04,path_MOYDIS04)
        ftp.quit()        
    MOD04 = res[0]    

    print 'load Aerosol band'
    if res=='no File':
	    tam=lat03_granule.shape
	    lat04_granule=lat03_granule[0:tam[0]:10,0:tam[1]:10]
	    lon04_granule=lon03_granule[0:tam[0]:10,0:tam[1]:10]
	    Aero= np.zeros(lon04_granule.shape)
	    Nl04, Ns04 =  lon04_granule.shape

    else:    
	    gdal_data = gdal.Open ('HDF4_EOS:EOS_SWATH:'+MOD04+':mod04:Deep_Blue_Aerosol_Optical_Depth_550_Land') 
	    gdal_data2 = gdal.Open ('HDF4_EOS:EOS_SWATH:'+MOD04+':mod04:AOD_550_Dark_Target_Deep_Blue_Combined') 
#	    gdal_data = gdal.Open ('HDF4_EOS:EOS_SWATH:'+MOD04+':mod04:Optical_Depth_Land_And_Ocean') 
	    Aero = np.array(gdal_data.ReadAsArray()[:,:],dtype=np.float) * float(gdal_data.GetMetadata()['scale_factor'].split(',')[0])
	    Aero2 = np.array(gdal_data2.ReadAsArray()[:,:],dtype=np.float) * float(gdal_data2.GetMetadata()['scale_factor'].split(',')[0])
            Aero[Aero<0]=Aero2[Aero<0]

	    gdal_data3 = gdal.Open ('HDF4_EOS:EOS_SWATH:'+MOD04+':mod04:Deep_Blue_Aerosol_Optical_Depth_550_Land_Estimated_Uncertainty')

	    AeroUn = np.array(gdal_data3.ReadAsArray()[:,:],dtype=np.float) * float(gdal_data3.GetMetadata()['scale_factor'].split(',')[0])
            #print AeroUn.shape, sum(AeroUn)

     	    gdal_data_lat = gdal.Open ('HDF4_EOS:EOS_SWATH:'+MOD04+':mod04:Latitude')
            lat04_granule = gdal_data_lat.ReadAsArray()
            gdal_data_lon = gdal.Open ('HDF4_EOS:EOS_SWATH:'+MOD04+':mod04:Longitude')
            lon04_granule = gdal_data_lon.ReadAsArray()
            Nl04, Ns04 =  lon04_granule.shape # Nbre line and sample
#	    print Nl04, Ns04


    if AOT_OMI_tag==1:



    ######################################## 
    #upload AURA OMI AOT product
    ######################################## 
	    OMAERO_tmp_T = 'OMI-Aura_L2G-OMAEROG_' + str(year) + 'm' + str(month).zfill(2) + str(day).zfill(2) +'_v003-*.he5'
	    res = glob.glob( path_OMAERO + OMAERO_tmp_T)
	    if len(res) == 0 :
		ftp = ftplib.FTP("acdisc.gsfc.nasa.gov")
		ftp.login("anonymous", "bernardo.mota@kcl.ac.uk")
		res = mt.download_atAura(ftp, OMAERO_tmp_T, year, doy, dir_OMAERO, path_OMAERO)
		ftp.quit()        
	    OMAE = res[0]    

	    print 'load OMI Aerosol data'
	#    if res=='no File':
	#    	    Aero= np.zeros((720,1440))#, dtype=np.uint8)
	#	    lat05_granule=lat03_granule[0:tam[0]:10,0:tam[1]:10]
	#	    lon05_granule=lon03_granule[0:tam[0]:10,0:tam[1]:10]
	#	    Aero= np.zeros(lon04_granule.shape)
	#	    Nl04, Ns04 =  lon04_granule.shape
	#    else:    
	    if res != 'no File':
		    gdal_data_AOT2 = gdal.Open ('HDF5:'+OMAE+'://HDFEOS/SWATHS/ColumnAmountAerosol/Data_Fields/AerosolOpticalThicknessMW') 
		    Aero_TEMP = np.array(gdal_data_AOT2.ReadAsArray()[:,:,:,:],dtype=np.float) * float(0.001)
		    Aero_OMI = Aero_TEMP[1,5,:,:]
	     	    gdal_data_lat = gdal.Open ('HDF5:'+OMAE+'://HDFEOS/SWATHS/ColumnAmountAerosol/Geolocation_Fields/Latitude')
		    lat05_granule = np.array(gdal_data_lat.ReadAsArray()[:,:],dtype=np.float) * float(1.0)
		    gdal_data_lon = gdal.Open ('HDF5:'+OMAE+'://HDFEOS/SWATHS/ColumnAmountAerosol/Geolocation_Fields/Longitude') 
		    lon05_granule = np.array(gdal_data_lon.ReadAsArray()[:,:],dtype=np.float) * float(1.0)
		    Nl05, Ns05 =  lon05_granule.shape # Nbre line and sample

    else:
	   Aero_OMI = Aero
	   lat05_granule = lat04_granule
           lon05_granule = lon04_granule
           Nl05, Ns05 =  Nl04, Ns04 


    ################################
    #compute fire characteristic
    ################################
    out, trans = mt.compute_cluster_info(fire_name,fire_product,fire_lat,fire_lon,zenith,azimuth,lat03_granule,lon03_granule,date_fire,date_modis,sat_,flag_ecmwf,ecmwf_dir,obs_hbl)

    #apply transmittance
    L_mir = L_mir/trans[0,0] 
    L_tir = L_tir/trans[1,0]

    # <headingcell level=1>
    ############################################################################################
    # 2.select a Region Of Interest
    ############################################################################################
    # <rawcell>
    # create ROI of "N_pixel" pixels around the fire location, and convert ROI's pixels (lat,lon) extracted from MOD03 to UTM
    # <codecell>

    #conversion to utm
    ##################
    #for MOD03 latlon coordinate (500m)
    Nbre_pt = Nl03* Ns03
    e03_granule = np.zeros(Nbre_pt,dtype=np.double)
    n03_granule = np.zeros(Nbre_pt,dtype=np.double)

    #for MOD04 latlon coordinate (10km)
    Nbre_pt4 = Nl04* Ns04
    e04_granule = np.zeros(Nbre_pt4,dtype=np.double)
    n04_granule = np.zeros(Nbre_pt4,dtype=np.double)

    #for OMAERO latlon coordinate (~12.5km)
    Nbre_pt5 = Nl05* Ns05
    e05_granule = np.zeros(Nbre_pt5,dtype=np.double)
    n05_granule = np.zeros(Nbre_pt5,dtype=np.double)

    #print 'number of pixels', Nbre_pt, Nbre_pt4, Nbre_pt5

    #--------------------------------------------------------------------------------------------------------
    #print 'convertion to utm'
    #fire location in the granule 
    dist_tmp = np.sqrt( (lat03_granule - fire_lat)**2 +  (lon03_granule - fire_lon)**2)
    #dist_tmp = np.sqrt( (lat03_granule - fire_lat[0])**2 +  (lon03_granule - fire_lon[0])**2)
    r_granule,c_granule = np.where( dist_tmp == dist_tmp.min() ) #Determine the centre of granule
    range_ii = np.arange( max([r_granule[0]-N_pixel,0]), min([r_granule[0]+N_pixel,Nl03]) )
    range_jj = np.arange( max([c_granule[0]-N_pixel,0]), min([c_granule[0]+N_pixel,Ns03]) )

    en_2_ijlatlon = np.zeros([Nbre_pt,2], dtype = np.int) # lookup table link between eastnorthing and latlon
    kk = 0
    for ii in range_ii:
        for jj in range_jj:
            e, n, z = conv_ll2utm.TransformPoint(float(lon03_granule[ii,jj]), float(lat03_granule[ii,jj]))
            e03_granule[kk],n03_granule[kk] = [e,n]
            en_2_ijlatlon[kk,:] = [ii,jj]
            kk += 1

    #--------------------------------------------------------------------------------------------------------
    #print 'convertion to utm'
    #fire location in the granule for 10km for AOT MODIS 
    dist_tmp4 = np.sqrt( (lat04_granule - fire_lat)**2 +  (lon04_granule - fire_lon)**2)
    r_granule4,c_granule4 = np.where( dist_tmp4 == dist_tmp4.min() )

    #print 'center pixels', r_granule4,c_granule4

    range_ii4 = np.arange( max([r_granule4[0]-N_pixel,0]), min([r_granule4[0]+N_pixel,Nl04]) )
    range_jj4 = np.arange( max([c_granule4[0]-N_pixel,0]), min([c_granule4[0]+N_pixel,Ns04]) )

    #print 'lines, colums,: 'range_ii4.shape, range_ii4.shape
    #print 'colunas: ', range_ii4.shape

    en_2_ijlatlon4 = np.zeros([Nbre_pt4,2], dtype = np.int) # lookup table link between eastnorthing and latlon
    kk = 0
    for ii in range_ii4:
        for jj in range_jj4:
            e, n, z = conv_ll2utm.TransformPoint(float(lon04_granule[ii,jj]), float(lat04_granule[ii,jj]))
            e04_granule[kk],n04_granule[kk] = [e,n]
            en_2_ijlatlon4[kk,:] = [ii,jj]
            kk += 1

    #---------------------------------------------------------------------------------------------------------
    #print 'convertion to utm'
    #fire location in the granule for 12.5km for AOT OMI resolutions
    dist_tmp5 = np.sqrt( (lat05_granule - fire_lat)**2 +  (lon05_granule - fire_lon)**2)
    r_granule5,c_granule5 = np.where( dist_tmp5 == dist_tmp5.min() )
    range_ii5 = np.arange( max([r_granule5[0]-N_pixel,0]), min([r_granule5[0]+N_pixel,Nl05]) )
    range_jj5 = np.arange( max([c_granule5[0]-N_pixel,0]), min([c_granule5[0]+N_pixel,Ns05]) )

    en_2_ijlatlon5 = np.zeros([Nbre_pt5,2], dtype = np.int) # lookup table link between eastnorthing and latlon
    kk = 0
    for ii in range_ii5:
        for jj in range_jj5:
            e, n, z = conv_ll2utm.TransformPoint(float(lon05_granule[ii,jj]), float(lat05_granule[ii,jj]))
            e05_granule[kk],n05_granule[kk] = [e,n]
            en_2_ijlatlon5[kk,:] = [ii,jj]
            kk += 1

    #---------------------------------------------------------------------------------------------------------            
    # convertion (lon_fire,lat_fire) to utm
    e_fire, n_fire, z = conv_ll2utm.TransformPoint(fire_lon, fire_lat)
            
    # define homogeneous grid of the ROI where the bands are going to be interpolated
    grid_e, grid_n = mt.defineGrid(e_fire, n_fire, box_size=georef_roi_box_size,res=dxy) #box_size in km


    #only keep fire product pixel located in the ROI
    idx = []
    for i_pixel in range(fire_product.shape[0]):
        if (fire_product.samp[i_pixel] >= range_jj.min()) & (fire_product.samp[i_pixel] <= range_jj.max()) & \
           (fire_product.line[i_pixel] >= range_ii.min()) & (fire_product.line[i_pixel] <= range_ii.max())   :
            idx.append(i_pixel)

    fire_product_roi = fire_product[idx]
    fire_product_roi  =fire_product_roi.view(np.recarray)
    print len(idx), ' fire pixels are in the ROI' 

#    llcrnrlon, llcrnrlat, tmp = conv_utm2ll.TransformPoint(grid_e[0,0]-dxy/2, grid_n[0,0]-dxy/2)
#    urcrnrlon, urcrnrlat, tmp = conv_utm2ll.TransformPoint(grid_e[-1,-1]+dxy/2, grid_n[-1,-1]+dxy/2)

    fire_msg_product_roi  =fire_MSG.view(np.recarray)



    print 'start georef ...',

    ################################
    # Georef individual Image types
    ################################
    # keep point in the box
    index_kept = np.where( (e03_granule >=  grid_e.min()) & (e03_granule <=  grid_e.max()) &\
                           (n03_granule >=  grid_n.min()) & (n03_granule <=  grid_n.max())  )[0]
    index_kept4 = np.where( (e04_granule >=  grid_e.min()) & (e04_granule <=  grid_e.max()) &\
                           (n04_granule >=  grid_n.min()) & (n04_granule <=  grid_n.max())  )[0]
    index_kept5 = np.where( (e05_granule >=  grid_e.min()) & (e05_granule <=  grid_e.max()) &\
                           (n05_granule >=  grid_n.min()) & (n05_granule <=  grid_n.max())  )[0]

#    print 'number of cases', index_kept, index_kept4, index_kept5

    # define variable
    Nbre_pt_kept  = len(index_kept) 
    data_easting  = np.zeros(Nbre_pt_kept)
    data_northing = np.zeros(Nbre_pt_kept)
    data          = np.zeros([Nbre_pt_kept,5])

    Nbre_pt_kept4  = len(index_kept4) 
    data_easting4  = np.zeros(Nbre_pt_kept4)
    data_northing4 = np.zeros(Nbre_pt_kept4)
    data4          = np.zeros([Nbre_pt_kept4,1])
    data4Un          = np.zeros([Nbre_pt_kept4,1])

    Nbre_pt_kept5  = len(index_kept5) 
    data_easting5  = np.zeros(Nbre_pt_kept5)
    data_northing5 = np.zeros(Nbre_pt_kept5)
    data5          = np.zeros([Nbre_pt_kept5,1])


#    print 'set data for lut.'
    if Nbre_pt_kept>0:
	percent_obs[0]=(Nbre_pt_kept*1.0)/(windowsize*windowsize/4.0)*100.0

    if Nbre_pt_kept4>0:
	percent_obs[1]=(Nbre_pt_kept4*1.0)/(windowsize*windowsize/400.0)*100.0

    if Nbre_pt_kept5>0:
	percent_obs[2]=(Nbre_pt_kept5*1.0)/(windowsize*windowsize/480.0)*100.0

    if Nbre_pt_kept == 0: 
        return out, grid_e, grid_n, dxy, conv_ll2utm, conv_utm2ll, AIcha, AIUncha, AISta, AIUnSta, Aero_geo, percent_obs
    #print 'set data for lut.', Nbre_pt_kept4
    if Nbre_pt_kept4 == 0: 
        return out, grid_e, grid_n, dxy, conv_ll2utm, conv_utm2ll, AIcha, AIUncha, AISta, AIUnSta, Aero_geo, percent_obs
    #print 'set data for lut.', Nbre_pt_kept5
    if Nbre_pt_kept5 == 0: 
        return out, grid_e, grid_n, dxy, conv_ll2utm, conv_utm2ll, AIcha, AIUncha, AISta, AIUnSta, Aero_geo, percent_obs

    data_easting   = e03_granule[index_kept]
    data_northing  = n03_granule[index_kept]
    data_easting4   = e04_granule[index_kept4]
    data_northing4  = n04_granule[index_kept4]
    data_easting5   = e05_granule[index_kept5]
    data_northing5  = n05_granule[index_kept5]

#    print 'make lut.'
    # georef False Color Composite
    for kk in range(Nbre_pt_kept):
        ii, jj  = en_2_ijlatlon[index_kept[kk],:]
        data[kk,0]     = L_mir[ii,jj]
        data[kk,1]     = L_tir[ii,jj]
        data[kk,2]     = imColorComp[ii,jj,0]
        data[kk,3]     = imColorComp[ii,jj,1]
        data[kk,4]     = imColorComp[ii,jj,2]

    for kk in range(Nbre_pt_kept4):
        ii, jj  = en_2_ijlatlon4[index_kept4[kk],:]
        data4[kk,0]     = Aero[ii,jj]
        data4Un[kk,0]   = AeroUn[ii,jj]

    for kk in range(Nbre_pt_kept5):
        ii, jj  = en_2_ijlatlon5[index_kept5[kk],:]
        data5[kk,0]     = Aero_OMI[ii,jj]


#    print 'stack value lut.'
    coord_pts = np.dstack((data_easting,data_northing))
    coord_pts4 = np.dstack((data_easting4,data_northing4))
    coord_pts5 = np.dstack((data_easting5,data_northing5))

    print 'AOT processing'
    ##########
    # AOT
    ##########
    Aero_geo = interpolate.griddata(coord_pts4.reshape(-1,2) , data4[:,0], (grid_e, grid_n), method='nearest' )
    Aero_geoUn = interpolate.griddata(coord_pts4.reshape(-1,2) , data4Un[:,0], (grid_e, grid_n), method='nearest' )


    AIcha, AIUncha, AISta, AIUnSta = aot.cal_stat_AOT(Aero_geo, Aero_geoUn)
    print 'AOT finnish'


#    numfilesUn=len(Aero_geoUn[Aero_geoUn>0])
#    numfiles=len(Aero_geo[Aero_geo>0])
#    Sum_Aero_geo=Aero_geo[Aero_geo>0].sum()
#    Sum_Aero_geoUn=Aero_geoUn[Aero_geoUn>0].sum()

#    AOT_tab=np.zeros(4)
#    AOT_tab=np.vstack((numfiles, numfilesUn, Sum_Aero_geo, Sum_Aero_geoUn))

    Aero_OMI_geo = interpolate.griddata(coord_pts5.reshape(-1,2) , data5[:,0], (grid_e, grid_n),fill_value=0, method='nearest' )

    ##########
    # mir
    ##########
    #use nearest neighbourgh to conserve radiance
    L_mir_geo = interpolate.griddata(coord_pts.reshape(-1,2) , data[:,0], (grid_e, grid_n),fill_value=0, method='nearest' )
    #convert to brightness temperature
    wavelength = np.zeros(L_mir_geo.shape) + 3.959
    BT_mir_geo = mt.planck_temperature(wavelength,L_mir_geo)
    # and apply the inverse of a model spectral function directly on BT to get central BT
    BT_mir_geo = mt.apply_inv_SRF_model_on_BT(BT_mir_geo,'MIR')

    ##########
    # tir
    ##########
    #use nearest neighbourgh to conserve radiance
    L_tir_geo = interpolate.griddata(coord_pts.reshape(-1,2) , data[:,1], (grid_e, grid_n),fill_value=0, method='nearest' )
    #convert to brightness temperature
    wavelength = np.zeros(L_tir_geo.shape) + 11.009
    BT_tir_geo = mt.planck_temperature(wavelength,L_tir_geo)
    # and apply a inverse of a model spectral function directly on BT to get central BT
    BT_tir_geo = mt.apply_inv_SRF_model_on_BT(BT_tir_geo,'TIR')

    ##########
    # true composite image and select color band
    ##########
    imColorComp_geo = np.zeros([grid_e.shape[0],grid_e.shape[1],3] , dtype=np.float)
    imColorComp_geo[:,:,0] = interpolate.griddata(coord_pts.reshape(-1,2) , data[:,2], (grid_e, grid_n), fill_value=0, method='linear' )
    imColorComp_geo[:,:,1] = interpolate.griddata(coord_pts.reshape(-1,2) , data[:,3], (grid_e, grid_n), fill_value=0, method='linear' )
    imColorComp_geo[:,:,2] = interpolate.griddata(coord_pts.reshape(-1,2) , data[:,4], (grid_e, grid_n), fill_value=0, method='linear' )    
    R = imColorComp_geo[:,:,0]; G = imColorComp_geo[:,:,1]; B = imColorComp_geo[:,:,2] # select color band
    imColorComp_geo_out = np.zeros(imColorComp_geo.shape)
    imColorComp_geo_out[:,:,0] = img_scale.sqrt(R, scale_min=0.025, scale_max=0.603)
    imColorComp_geo_out[:,:,1] = img_scale.sqrt(G, scale_min=0.044, scale_max=0.546)
    imColorComp_geo_out[:,:,2] = img_scale.sqrt(B, scale_min=0.082, scale_max=0.525)

    print 'done'

    if tagfigures==1:
	#############################################################################################################################
	# print MIR-TIR raw difference imagery
	#############################################################################################################################

	    if tagprint1==1:

		    # <headingcell level=1>
		    ############################################################################################
		    # 5. Plot georeference Data
		    ############################################################################################
		    
		    # <headingcell level=2>
		    # 5.1 plot MIR and TIR
		    ############################################################################################
		    # <codecell>
		    mpl.rcdefaults()
		    fig = plt.figure(figsize=(20,40))
		    extent = (grid_e.min()-dxy/2,grid_e.max()+dxy/2,grid_n.min()-dxy/2,grid_n.max()+dxy/2)

		    ax = plt.subplot(311)
		    im = ax.imshow(BT_mir_geo.T,origin='lower',interpolation='nearest',extent=extent)
		    cbar = plt.colorbar(im)
		    cbar.set_label('BT (K)')
		    ax.set_axis_off()
		    ax.set_title('BT MIR')

		    ax = plt.subplot(312)
		    im = ax.imshow(BT_tir_geo.T,origin='lower',interpolation='nearest',extent=extent)
		    cbar = plt.colorbar(im)
		    cbar.set_label('BT (K)')
		    ax.set_axis_off()
		    ax.set_title('BT TIR')

		    ax = plt.subplot(313)
		    diff = L_mir_geo - L_tir_geo
		    im = ax.imshow(diff.T,origin='lower',interpolation='nearest',extent=extent)
		    cbar = plt.colorbar(im)
		    cbar.set_label('Delta radiance (W/m2/str)')
		    ax.set_axis_off()
		    ax.set_title('Radiance: MIR - TIR')

		    #plot hot spot
		    idx = []
		    for i_pixel in range(fire_product_roi.shape[0]):
			if (fire_product_roi.northing[i_pixel] >= grid_n.min()) & (fire_product_roi.northing[i_pixel] <= grid_n.max()) & \
			   (fire_product_roi.easting[i_pixel] >= grid_e.min()) & (fire_product_roi.easting[i_pixel] <= grid_e.max())   :
			    idx.append(i_pixel)

		    ax.scatter(fire_product_roi.easting[idx],fire_product_roi.northing[idx],\
			       edgecolors='k',marker='+',facecolors='None',s=70,alpha=.7,linewidth='3')
		    ax.set_xlim(extent[0],extent[1])
		    ax.set_ylim(extent[2],extent[3])

		    fig.savefig(path_image+'geoRef_MIR_TIR_'+'{:04d}{:03d}{:04d}_{:s}.png'.format(year,doy,int(float(time)),sat))
		    plt.close(fig)


	#############################################################################################################################
	# print geoRef_MIR image
	#############################################################################################################################


	    #plot only MIR
	    mpl.rcdefaults()
	    mpl.rcParams['text.usetex'] = True
	    mpl.rcParams['font.family'] = 'Comic Sans MS'
	    mpl.rcParams['axes.linewidth'] = 1
	    mpl.rcParams['axes.labelsize'] = 24.
	    mpl.rcParams['legend.fontsize'] = 'small'
	    mpl.rcParams['legend.fancybox'] = True
	    mpl.rcParams['font.size'] = 24.
	    mpl.rcParams['xtick.labelsize'] = 24.
	    mpl.rcParams['ytick.labelsize'] = 24.
	    mpl.rcParams['figure.subplot.left'] = .0
	    mpl.rcParams['figure.subplot.right'] = 1.
	    mpl.rcParams['figure.subplot.top'] = 1.
	    mpl.rcParams['figure.subplot.bottom'] = .0
	    mpl.rcParams['figure.subplot.hspace'] = 0.1
	    mpl.rcParams['figure.subplot.wspace'] = 0.18
	    fig = plt.figure(figsize=(12,12))
	    extent = (grid_e.min()-dxy/2,grid_e.max()+dxy/2,grid_n.min()-dxy/2,grid_n.max()+dxy/2)

	    ax = plt.subplot(111)
	    
	    # setup of basemap ('lcc' = lambert conformal conic).
	    # use major and minor sphere radii from WGS84 ellipsoid.
	    nx, ny = grid_e.shape
	    llcrnrlon, llcrnrlat, tmp = conv_utm2ll.TransformPoint(grid_e[0,0]      -dxy/2, grid_n[0,0]     -dxy/2)
	    urcrnrlon, urcrnrlat, tmp = conv_utm2ll.TransformPoint(grid_e[-1,-1]    +dxy/2, grid_n[-1,-1]   +dxy/2)
	    width_ =   (grid_e[-1,-1]    +dxy/2) - (grid_e[0,0]      -dxy/2)
	    height_ =  (grid_n[-1,-1]    +dxy/2) - (grid_n[0,0]      -dxy/2)

	    lon_0, lat_0, tmp         = conv_utm2ll.TransformPoint(grid_e[nx/2,ny/2],grid_n[nx/2,ny/2])
	    #m = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, projection='tmerc', lon_0=lon_0, lat_0=lat_0, resolution='i')
	    m = Basemap(width=width_, height=height_, projection='tmerc', lon_0=lon_0, lat_0=lat_0, resolution='i')

	    extent=(llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat)
	    # plot image over map with imshow.
	    #cmap_h = cmplt.get_cmap('Dark2',9 )
	    
	    #m.drawcoastlines()
	    #m.fillcontinents(color='coral',lake_color='aqua')
	    # draw parallels and meridians.
	    #m.drawparallels(np.arange(-40,61.,2.))
	    #m.drawmeridians(np.arange(-20.,21.,2.))
	    #m.drawmapboundary(fill_color='aqua')
	    im = m.imshow(np.transpose(BT_mir_geo),interpolation='nearest',extent=extent,cmap=mpl.cm.CMRmap)
	    #plt.show()
	    #pdb.set_trace()

	    # draw coastlines and political boundaries.
	    #m.drawcoastlines()
	    #m.drawcountries()
	    #m.drawstates() 
	    
	    #plot hot spot
	    idx = []
	    for i_pixel in range(fire_product_roi.shape[0]):
		if (fire_product_roi.northing[i_pixel] >= grid_n.min()) & (fire_product_roi.northing[i_pixel] <= grid_n.max()) & \
		   (fire_product_roi.easting[i_pixel] >= grid_e.min()) & (fire_product_roi.easting[i_pixel] <= grid_e.max())   :
		    idx.append(i_pixel)

	    idx_msg = []
	    for i_pixel in range(fire_msg_product_roi.shape[0]):
		if (fire_msg_product_roi.northing[i_pixel] >= grid_n.min()) & (fire_msg_product_roi.northing[i_pixel] <= grid_n.max()) & \
		   (fire_msg_product_roi.easting[i_pixel] >= grid_e.min()) & (fire_msg_product_roi.easting[i_pixel] <= grid_e.max())   :
		    idx_msg.append(i_pixel)


	    # Draw a map scale
	    dev = m.scatter(fire_product_roi.lon[idx],fire_product_roi.lat[idx],latlon=True,
		marker='+', lw=3., s=100, 
		facecolor='g', edgecolor='g',
		alpha=1., antialiased=True,
		label='MOD14 hot spot', zorder=3)

	    dev = m.scatter(fire_msg_product_roi.lon[idx_msg],fire_msg_product_roi.lat[idx_msg],latlon=True,
		marker='+', lw=3., s=100, 
		facecolor='w', edgecolor='w',
		alpha=1., antialiased=True,
		label='MSG hot spot', zorder=3)

	    if idx_msg>0:
	    	total_msg_frp=sum(fire_msg_product_roi.FRP[idx_msg])
	    else:
		total_msg_frp=0.0

	    m.drawmapscale(
		llcrnrlon+0.15, llcrnrlat+.05 + 0.015,
		urcrnrlon-0.25,        llcrnrlat+.05 ,
		10.,
		barstyle='fancy', labelstyle='simple',
		fillcolor1='w', fillcolor2='k',
		fontcolor='w',fontsize=24,
		zorder=5)
	    
	    #add compass
	    north_indicator =  np.array(Image.open('../data_static/compass_invertedColor.png'))
	    im_compass = OffsetImage(north_indicator, zoom=.3)
	    ab = AnnotationBbox(im_compass, [.08*m.urcrnrx,.9*m.urcrnry], xycoords='data', frameon=False)
	    ax.add_artist(ab)
	    
	    #im = ax.imshow(BT_mir_geo.T,origin='lower',interpolation='nearest',extent=extent,cmap=mpl.cm.CMRmap)
	    #divider = make_axes_locatable(ax)
	    #cbaxes = divider.append_axes("right", size="5%", pad=0.05)
	    cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.8])
	    cbar = fig.colorbar(im,orientation='vertical',cax = cbaxes)
	    
	    cbytick_obj = plt.getp(cbar.ax.axes, 'yticklabels') 
	    plt.setp(cbytick_obj, color='w')
	    cbar.outline.set_color('w')                   #set colorbar box color
	    cbar.ax.yaxis.set_tick_params(color='w')      #set colorbar ticks color 
	    
	    #cbar = m.colorbar(im,location='right', size="5%", pad=-0.05)
	    #cbar = fig.colorbar(im ,cax = cbaxes)
	    cbar.set_label('Brightness Temperature (K)',labelpad=10,color='w')

	    ax.set_axis_off()
	    ax.set_title(r'{:04d}-{:02d}-{:02d} - {:02d}:{:02d} - {:s} vz={:3.1f} va={:3.1f}'.format(date_modis.year,date_modis.month,date_modis.day,\
		                                                                          date_modis.hour,date_modis.minute,sat_,out.viewAngle[0],out.azimuthAngle[0]),\
		                                                                          y=.96,color='w')

	    #add FRP info
	    # these are matplotlib.patch.Patch properties
	    
	    props = dict(boxstyle='round', facecolor='white', alpha=0.5)

	    # place a text box in upper left in axes coords
	    textstr=''
	    print out.FRP_kcl_bF_MW[0]
	    print total_msg_frp
	    if out.FRP_kcl_bnF_MW[0] > 0:
		textstr  = '$FRP=%.2f ~MW$ \n'%(out.FRP_kcl_bnF_MW[0])
		textstr += '$AF area=%.2f ~ha$'%(out.afArea_bnF_ha[0])
	#    else:
	#	textstr  = '$FRP=%.2f ~MW$ \n'%(0.0)
	    textstr2=''
	    if out.FRP_kcl_bF_MW[0] > 0:
	 #      textstr2  = '$FRP^{f}=%.2f ~MW$ \n'%(total_msg_frp)
		textstr2  = '$FRP^{f}=%.2f ~MW$ \n'%(out.FRP_kcl_bF_MW[0])
		#textstr2 += '$AF area^{f}=%.2f ~ha$'%(out.afArea_bF_ha[0])
		textstr2 += '$MSG FRP=%.2f ~MW$'%(total_msg_frp)
	    else:
		textstr2  = '$FRP=%.2f ~MW$ \n'%(0.0)
		textstr2 += '$MSG FRP=%.2f ~MW$'%(total_msg_frp)

	    print textstr
	    print textstr2
	    if textstr != '':
		ax.text(0.2, 0.12, textstr, transform=ax.transAxes, fontsize=24, verticalalignment='top', bbox=props)
	    if textstr2 != '':
		ax.text(0.5, 0.12, textstr2, transform=ax.transAxes, fontsize=24, verticalalignment='top', bbox=props)
	    
	    #ax.scatter(fire_product_roi.easting[idx],fire_product_roi.northing[idx],\
	    #           edgecolors='g',marker='+',facecolors='g',s=70,alpha=1.,linewidth='1')
	    #ax.set_xlim(extent[0],extent[1])
	    #ax.set_ylim(extent[2],extent[3])

	    fig.savefig(path_image+'geoRef_MIR_'+'{:04d}{:03d}{:04d}_{:s}.png'.format(year,doy,int(float(time)),sat))
	    plt.close(fig)
	    #for the kml file
	    out.path_geoRef_MIR[0] = path_image+'geoRef_MIR_'+'{:04d}{:03d}{:04d}_{:s}.png'.format(year,doy,int(float(time)),sat)


	#############################################################################################################################
	# print MODIS Aerosol image
	#############################################################################################################################


	    #plot only Aero
	    mpl.rcdefaults()
	    mpl.rcParams['text.usetex'] = True
	    mpl.rcParams['font.family'] = 'Comic Sans MS'
	    mpl.rcParams['axes.linewidth'] = 1
	    mpl.rcParams['axes.labelsize'] = 24.
	    mpl.rcParams['legend.fontsize'] = 'small'
	    mpl.rcParams['legend.fancybox'] = True
	    mpl.rcParams['font.size'] = 24.
	    mpl.rcParams['xtick.labelsize'] = 24.
	    mpl.rcParams['ytick.labelsize'] = 24.
	    mpl.rcParams['figure.subplot.left'] = .0
	    mpl.rcParams['figure.subplot.right'] = 1.
	    mpl.rcParams['figure.subplot.top'] = 1.
	    mpl.rcParams['figure.subplot.bottom'] = .0
	    mpl.rcParams['figure.subplot.hspace'] = 0.1
	    mpl.rcParams['figure.subplot.wspace'] = 0.18
	    fig = plt.figure(figsize=(12,12))
	    extent = (grid_e.min()-dxy/2,grid_e.max()+dxy/2,grid_n.min()-dxy/2,grid_n.max()+dxy/2)

	    ax = plt.subplot(111)
	    
	    # setup of basemap ('lcc' = lambert conformal conic).
	    # use major and minor sphere radii from WGS84 ellipsoid.
	    nx, ny = grid_e.shape
	    llcrnrlon, llcrnrlat, tmp = conv_utm2ll.TransformPoint(grid_e[0,0]      -dxy/2, grid_n[0,0]     -dxy/2)
	    urcrnrlon, urcrnrlat, tmp = conv_utm2ll.TransformPoint(grid_e[-1,-1]    +dxy/2, grid_n[-1,-1]   +dxy/2)
	    width_ =   (grid_e[-1,-1]    +dxy/2) - (grid_e[0,0]      -dxy/2)
	    height_ =  (grid_n[-1,-1]    +dxy/2) - (grid_n[0,0]      -dxy/2)

	    lon_0, lat_0, tmp         = conv_utm2ll.TransformPoint(grid_e[nx/2,ny/2],grid_n[nx/2,ny/2])
	    #m = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, projection='tmerc', lon_0=lon_0, lat_0=lat_0, resolution='i')
	    m = Basemap(width=width_, height=height_, projection='tmerc', lon_0=lon_0, lat_0=lat_0, resolution='i')

	    extent=(llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat)
	    # plot image over map with imshow.
	    #cmap_h = cmplt.get_cmap('Dark2',9 )
	    
	    #m.drawcoastlines()
	    #m.fillcontinents(color='coral',lake_color='aqua')
	    # draw parallels and meridians.
	    #m.drawparallels(np.arange(-40,61.,2.))
	    #m.drawmeridians(np.arange(-20.,21.,2.))
	    #m.drawmapboundary(fill_color='aqua')
	    im = m.imshow(np.transpose(img_scale.linear(Aero_geo, scale_min=-0.1, scale_max=5.0)),interpolation='nearest',extent=extent,cmap=mpl.cm.CMRmap)
	    #plt.show()
	    #pdb.set_trace()

	    # draw coastlines and political boundaries.
	    #m.drawcoastlines()
	    #m.drawcountries()
	    #m.drawstates() 
	    
	    #plot hot spot
	    idx = []
	    for i_pixel in range(fire_product_roi.shape[0]):
		if (fire_product_roi.northing[i_pixel] >= grid_n.min()) & (fire_product_roi.northing[i_pixel] <= grid_n.max()) & \
		   (fire_product_roi.easting[i_pixel] >= grid_e.min()) & (fire_product_roi.easting[i_pixel] <= grid_e.max())   :
		    idx.append(i_pixel)

	    idx_msg = []
	    for i_pixel in range(fire_msg_product_roi.shape[0]):
		if (fire_msg_product_roi.northing[i_pixel] >= grid_n.min()) & (fire_msg_product_roi.northing[i_pixel] <= grid_n.max()) & \
		   (fire_msg_product_roi.easting[i_pixel] >= grid_e.min()) & (fire_msg_product_roi.easting[i_pixel] <= grid_e.max())   :
		    idx_msg.append(i_pixel)




	    # Draw a map scale
	    dev = m.scatter(fire_product_roi.lon[idx],fire_product_roi.lat[idx],latlon=True,
		marker='+', lw=3., s=100, 
		facecolor='g', edgecolor='g',
		alpha=1., antialiased=True,
		label='MOD14 hot spot', zorder=3)

	    dev = m.scatter(fire_msg_product_roi.lon[idx_msg],fire_msg_product_roi.lat[idx_msg],latlon=True,
		marker='+', lw=3., s=100, 
		facecolor='w', edgecolor='w',
		alpha=1., antialiased=True,
		label='MSG hot spot', zorder=3)

	    m.drawmapscale(
		llcrnrlon+0.15, llcrnrlat+.05 + 0.015,
		urcrnrlon-0.25,        llcrnrlat+.05 ,
		10.,
		barstyle='fancy', labelstyle='simple',
		fillcolor1='w', fillcolor2='k',
		fontcolor='w',fontsize=24,
		zorder=5)
	    
	    #add compass
	    north_indicator =  np.array(Image.open('../data_static/compass_invertedColor.png'))
	    im_compass = OffsetImage(north_indicator, zoom=.3)
	    ab = AnnotationBbox(im_compass, [.08*m.urcrnrx,.9*m.urcrnry], xycoords='data', frameon=False)
	    ax.add_artist(ab)
	    
	    #im = ax.imshow(BT_mir_geo.T,origin='lower',interpolation='nearest',extent=extent,cmap=mpl.cm.CMRmap)
	    #divider = make_axes_locatable(ax)
	    #cbaxes = divider.append_axes("right", size="5%", pad=0.05)
	    cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.8])
	    cbar = fig.colorbar(im,orientation='vertical',cax = cbaxes)
	    
	    cbytick_obj = plt.getp(cbar.ax.axes, 'yticklabels') 
	    plt.setp(cbytick_obj, color='b')
	    cbar.outline.set_color('b')                   #set colorbar box color
	    cbar.ax.yaxis.set_tick_params(color='b')      #set colorbar ticks color 
	    
	    #cbar = m.colorbar(im,location='right', size="5%", pad=-0.05)
	    #cbar = fig.colorbar(im ,cax = cbaxes)
	    cbar.set_label('AOT',labelpad=10,color='b')

	    ax.set_axis_off()
	    ax.set_title(r'{:04d}-{:02d}-{:02d} - {:02d}:{:02d} - {:s} vz={:3.1f} va={:3.1f}'.format(date_modis.year,date_modis.month,date_modis.day,\
		                                                                          date_modis.hour,date_modis.minute,sat_,out.viewAngle[0],out.azimuthAngle[0]),\
		                                                                          y=.96,color='w')

	    #add FRP info
	    # these are matplotlib.patch.Patch properties
	    
	    props = dict(boxstyle='round', facecolor='white', alpha=0.5)

	    # place a text box in upper left in axes coords
	    textstr=''
	    if out.FRP_kcl_bnF_MW[0] > 0:
		textstr  = '$FRP=%.2f ~MW$ \n'%(out.FRP_kcl_bnF_MW[0])
		textstr += '$AF area=%.2f ~ha$'%(out.afArea_bnF_ha[0])

	    textstr2=''
	    if out.FRP_kcl_bF_MW[0] > 0:
		textstr2  = '$FRP^{f}=%.2f ~MW$ \n'%(out.FRP_kcl_bF_MW[0])
		#textstr2 += '$AF area^{f}=%.2f ~ha$'%(out.afArea_bF_ha[0])
		#textstr2 += '$AF area^{f}=%.2f ~ha$'%(out.afArea_bF_ha[0])
		textstr2 += '$MSG FRP=%.2f ~MW$'%(total_msg_frp)
	    else:
		textstr2  = '$FRP=%.2f ~MW$ \n'%(0.0)
		textstr2 += '$MSG FRP=%.2f ~MW$'%(total_msg_frp)

	    if textstr != '':
		ax.text(0.2, 0.12, textstr, transform=ax.transAxes, fontsize=24, verticalalignment='top', bbox=props)
	    if textstr2 != '':
		ax.text(0.5, 0.12, textstr2, transform=ax.transAxes, fontsize=24, verticalalignment='top', bbox=props)
	    
	    #ax.scatter(fire_product_roi.easting[idx],fire_product_roi.northing[idx],\
	    #           edgecolors='g',marker='+',facecolors='g',s=70,alpha=1.,linewidth='1')
	    #ax.set_xlim(extent[0],extent[1])
	    #ax.set_ylim(extent[2],extent[3])

	    fig.savefig(path_image+'geoRef_AOT_'+'{:04d}{:03d}{:04d}_{:s}.png'.format(year,doy,int(float(time)),sat))
	    plt.close(fig)
	    #for the kml file
	    #out.path_geoRef_AOT[0] = path_image+'geoRef_AOT_'+'{:04d}{:03d}{:04d}_{:s}.png'.format(year,doy,int(float(time)),sat)


	#############################################################################################################################
	# print AURA OMI Aerosol image
	#############################################################################################################################

	    if AOT_OMI_tag==1:
		    #plot only Aero
		    mpl.rcdefaults()
		    mpl.rcParams['text.usetex'] = True
		    mpl.rcParams['font.family'] = 'Comic Sans MS'
		    mpl.rcParams['axes.linewidth'] = 1
		    mpl.rcParams['axes.labelsize'] = 24.
		    mpl.rcParams['legend.fontsize'] = 'small'
		    mpl.rcParams['legend.fancybox'] = True
		    mpl.rcParams['font.size'] = 24.
		    mpl.rcParams['xtick.labelsize'] = 24.
		    mpl.rcParams['ytick.labelsize'] = 24.
		    mpl.rcParams['figure.subplot.left'] = .0
		    mpl.rcParams['figure.subplot.right'] = 1.
		    mpl.rcParams['figure.subplot.top'] = 1.
		    mpl.rcParams['figure.subplot.bottom'] = .0
		    mpl.rcParams['figure.subplot.hspace'] = 0.1
		    mpl.rcParams['figure.subplot.wspace'] = 0.18
		    fig = plt.figure(figsize=(12,12))
		    extent = (grid_e.min()-dxy/2,grid_e.max()+dxy/2,grid_n.min()-dxy/2,grid_n.max()+dxy/2)

		    ax = plt.subplot(111)
		    
		    # setup of basemap ('lcc' = lambert conformal conic).
		    # use major and minor sphere radii from WGS84 ellipsoid.
		    nx, ny = grid_e.shape
		    llcrnrlon, llcrnrlat, tmp = conv_utm2ll.TransformPoint(grid_e[0,0]      -dxy/2, grid_n[0,0]     -dxy/2)
		    urcrnrlon, urcrnrlat, tmp = conv_utm2ll.TransformPoint(grid_e[-1,-1]    +dxy/2, grid_n[-1,-1]   +dxy/2)
		    width_ =   (grid_e[-1,-1]    +dxy/2) - (grid_e[0,0]      -dxy/2)
		    height_ =  (grid_n[-1,-1]    +dxy/2) - (grid_n[0,0]      -dxy/2)

		    lon_0, lat_0, tmp         = conv_utm2ll.TransformPoint(grid_e[nx/2,ny/2],grid_n[nx/2,ny/2])
		    #m = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, projection='tmerc', lon_0=lon_0, lat_0=lat_0, resolution='i')
		    m = Basemap(width=width_, height=height_, projection='tmerc', lon_0=lon_0, lat_0=lat_0, resolution='i')

		    extent=(llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat)
		    # plot image over map with imshow.
		    #cmap_h = cmplt.get_cmap('Dark2',9 )
		    
		    #m.drawcoastlines()
		    #m.fillcontinents(color='coral',lake_color='aqua')
		    # draw parallels and meridians.
		    #m.drawparallels(np.arange(-40,61.,2.))
		    #m.drawmeridians(np.arange(-20.,21.,2.))
		    #m.drawmapboundary(fill_color='aqua')
		    im = m.imshow(np.transpose(img_scale.sqrt(Aero_OMI_geo, scale_min=-0.1, scale_max=5.0)),interpolation='nearest',extent=extent,cmap=mpl.cm.CMRmap)
		    #plt.show()
		    #pdb.set_trace()

		    # draw coastlines and political boundaries.
		    #m.drawcoastlines()
		    #m.drawcountries()
		    #m.drawstates() 
		    
		    #plot hot spot
		    idx = []
		    for i_pixel in range(fire_product_roi.shape[0]):
			if (fire_product_roi.northing[i_pixel] >= grid_n.min()) & (fire_product_roi.northing[i_pixel] <= grid_n.max()) & \
			   (fire_product_roi.easting[i_pixel] >= grid_e.min()) & (fire_product_roi.easting[i_pixel] <= grid_e.max())   :
			    idx.append(i_pixel)

		    idx_msg = []
		    for i_pixel in range(fire_msg_product_roi.shape[0]):
			if (fire_msg_product_roi.northing[i_pixel] >= grid_n.min()) & (fire_msg_product_roi.northing[i_pixel] <= grid_n.max()) & \
			   (fire_msg_product_roi.easting[i_pixel] >= grid_e.min()) & (fire_msg_product_roi.easting[i_pixel] <= grid_e.max())   :
			    idx_msg.append(i_pixel)




		    # Draw a map scale
		    dev = m.scatter(fire_product_roi.lon[idx],fire_product_roi.lat[idx],latlon=True,
			marker='+', lw=3., s=100, 
			facecolor='g', edgecolor='g',
			alpha=1., antialiased=True,
			label='MOD14 hot spot', zorder=3)

		    dev = m.scatter(fire_msg_product_roi.lon[idx_msg],fire_msg_product_roi.lat[idx_msg],latlon=True,
			marker='+', lw=3., s=100, 
			facecolor='w', edgecolor='w',
			alpha=1., antialiased=True,
			label='MSG hot spot', zorder=3)

		    m.drawmapscale(
			llcrnrlon+0.15, llcrnrlat+.05 + 0.015,
			urcrnrlon-0.25,        llcrnrlat+.05 ,
			10.,
			barstyle='fancy', labelstyle='simple',
			fillcolor1='w', fillcolor2='k',
			fontcolor='w',fontsize=24,
			zorder=5)
		    
		    #add compass
		    north_indicator =  np.array(Image.open('../data_static/compass_invertedColor.png'))
		    im_compass = OffsetImage(north_indicator, zoom=.3)
		    ab = AnnotationBbox(im_compass, [.08*m.urcrnrx,.9*m.urcrnry], xycoords='data', frameon=False)
		    ax.add_artist(ab)
		    
		    #im = ax.imshow(BT_mir_geo.T,origin='lower',interpolation='nearest',extent=extent,cmap=mpl.cm.CMRmap)
		    #divider = make_axes_locatable(ax)
		    #cbaxes = divider.append_axes("right", size="5%", pad=0.05)
		    cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.8])
		    cbar = fig.colorbar(im,orientation='vertical',cax = cbaxes)
		    
		    cbytick_obj = plt.getp(cbar.ax.axes, 'yticklabels') 
		    plt.setp(cbytick_obj, color='b')
		    cbar.outline.set_color('b')                   #set colorbar box color
		    cbar.ax.yaxis.set_tick_params(color='b')      #set colorbar ticks color 
		    
		    #cbar = m.colorbar(im,location='right', size="5%", pad=-0.05)
		    #cbar = fig.colorbar(im ,cax = cbaxes)
		    cbar.set_label('AOT',labelpad=10,color='b')

		    ax.set_axis_off()
		    ax.set_title(r'{:04d}-{:02d}-{:02d} - {:02d}:{:02d} - {:s} vz={:3.1f} va={:3.1f}'.format(date_modis.year,date_modis.month,date_modis.day,\
				                                                                  date_modis.hour,date_modis.minute,sat_,out.viewAngle[0],out.azimuthAngle[0]),\
				                                                                  y=.96,color='w')

		    #add FRP info
		    # these are matplotlib.patch.Patch properties
		    
		    props = dict(boxstyle='round', facecolor='white', alpha=0.5)

		    # place a text box in upper left in axes coords
		    textstr=''
		    if out.FRP_kcl_bnF_MW[0] > 0:
			textstr  = '$FRP=%.2f ~MW$ \n'%(out.FRP_kcl_bnF_MW[0])
			textstr += '$AF area=%.2f ~ha$'%(out.afArea_bnF_ha[0])

		    textstr2=''
		    if out.FRP_kcl_bF_MW[0] > 0:
			textstr2  = '$FRP^{f}=%.2f ~MW$ \n'%(out.FRP_kcl_bF_MW[0])
			#textstr2 += '$AF area^{f}=%.2f ~ha$'%(out.afArea_bF_ha[0])
			#textstr2 += '$AF area^{f}=%.2f ~ha$'%(out.afArea_bF_ha[0])
			textstr2 += '$MSG FRP=%.2f ~MW$'%(total_msg_frp)
		    else:
			textstr2  = '$FRP=%.2f ~MW$ \n'%(0.0)
			textstr2 += '$MSG FRP=%.2f ~MW$'%(total_msg_frp)

		    if textstr != '':
			ax.text(0.2, 0.12, textstr, transform=ax.transAxes, fontsize=24, verticalalignment='top', bbox=props)
		    if textstr2 != '':
			ax.text(0.5, 0.12, textstr2, transform=ax.transAxes, fontsize=24, verticalalignment='top', bbox=props)
		    
		    #ax.scatter(fire_product_roi.easting[idx],fire_product_roi.northing[idx],\
		    #           edgecolors='g',marker='+',facecolors='g',s=70,alpha=1.,linewidth='1')
		    #ax.set_xlim(extent[0],extent[1])
		    #ax.set_ylim(extent[2],extent[3])

		    fig.savefig(path_image+'geoRef_OMI_AOT_'+'{:04d}{:03d}{:04d}_{:s}.png'.format(year,doy,int(float(time)),sat))
		    plt.close(fig)
		    #for the kml file
		    #out.path_geoRef_AOT[0] = path_image+'geoRef_AOT_'+'{:04d}{:03d}{:04d}_{:s}.png'.format(year,doy,int(float(time)),sat)



	#############################################################################################################################
	# print MODIS Aerosol image
	#############################################################################################################################


	    # <headingcell level=2>
	    # 5.2 Plot True Color Composite Image
	    ############################################################################################
	    # <codecell>

	    mpl.rcdefaults()
	    mpl.rcParams['text.usetex'] = True
	    mpl.rcParams['font.family'] = 'Comic Sans MS'
	    mpl.rcParams['axes.linewidth'] = 1
	    mpl.rcParams['axes.labelsize'] = 24.
	    mpl.rcParams['legend.fontsize'] = 'small'
	    mpl.rcParams['legend.fancybox'] = True
	    mpl.rcParams['font.size'] = 24.
	    mpl.rcParams['xtick.labelsize'] = 24.
	    mpl.rcParams['ytick.labelsize'] = 24.
	    mpl.rcParams['figure.subplot.left'] = .0
	    mpl.rcParams['figure.subplot.right'] = 1.
	    mpl.rcParams['figure.subplot.top'] = 1.
	    mpl.rcParams['figure.subplot.bottom'] = .0
	    mpl.rcParams['figure.subplot.hspace'] = 0.1
	    mpl.rcParams['figure.subplot.wspace'] = 0.18
	    fig = plt.figure(figsize=(12,12))
	    ax = plt.subplot(111)
	    # setup of basemap ('lcc' = lambert conformal conic).
	    # use major and minor sphere radii from WGS84 ellipsoid.
	    nx, ny = grid_e.shape
	    llcrnrlon, llcrnrlat, tmp = conv_utm2ll.TransformPoint(grid_e[0,0]      -dxy/2, grid_n[0,0]     -dxy/2)
	    urcrnrlon, urcrnrlat, tmp = conv_utm2ll.TransformPoint(grid_e[-1,-1]    +dxy/2, grid_n[-1,-1]   +dxy/2)
	    width_ =   (grid_e[-1,-1]    +dxy/2) - (grid_e[0,0]      -dxy/2)
	    height_ =  (grid_n[-1,-1]    +dxy/2) - (grid_n[0,0]      -dxy/2)

	    lon_0, lat_0, tmp         = conv_utm2ll.TransformPoint(grid_e[nx/2,ny/2],grid_n[nx/2,ny/2])
	    #m = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, projection='tmerc', lon_0=lon_0, lat_0=lat_0, resolution='i')
	    m = Basemap(width=width_, height=height_, projection='tmerc', lon_0=lon_0, lat_0=lat_0, resolution='i')

	    extent=(llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat)
	    # plot image over map with imshow.
	    #cmap_h = cmplt.get_cmap('Dark2',9 )
	    im = m.imshow(np.transpose(imColorComp_geo_out,[1,0,2]),interpolation='nearest',extent=extent)

	    # draw coastlines and political boundaries.
	    #m.drawcoastlines()
	    #m.drawcountries()
	    #m.drawstates() 

	    # Draw a map scale
	    dev = m.scatter(fire_product_roi.lon[idx],fire_product_roi.lat[idx],latlon=True,
		marker='+', lw=2., s=100,
		facecolor='r', edgecolor='r',
		alpha=1., antialiased=True,
		label='MOD14 hot spot', zorder=3)

	    dev = m.scatter(fire_msg_product_roi.lon[idx_msg],fire_msg_product_roi.lat[idx_msg],latlon=True,
		marker='+', lw=3., s=100, 
		facecolor='w', edgecolor='w',
		alpha=1., antialiased=True,
		label='MSG hot spot', zorder=3)

	    m.drawmapscale(
		urcrnrlon-0.25 + 0.08, llcrnrlat+.05 + 0.015,
		urcrnrlon-0.25,        llcrnrlat+.05 ,
		10.,
		barstyle='fancy', labelstyle='simple',
		fillcolor1='w', fillcolor2='k',
		fontcolor='k',fontsize=24,
		zorder=5)
	    
	    #add compass
	    north_indicator =  np.array(Image.open('../data_static/compass_invertedColor.png'))
	    im = OffsetImage(north_indicator, zoom=.3)
	    ab = AnnotationBbox(im, [.08*m.urcrnrx,.9*m.urcrnry], xycoords='data', frameon=False)
	    ax.add_artist(ab)

	    #add time
	    ax.set_title(r'{:04d}-{:02d}-{:02d} - {:02d}:{:02d} - {:s} vz={:3.1f} va={:3.1f}'.\
		         format(date_modis.year,date_modis.month,date_modis.day,date_modis.hour,date_modis.minute,sat_,out.viewAngle[0],out.azimuthAngle[0]), y=.96)

	    #add cluster info
	    # these are matplotlib.patch.Patch properties
	    props = dict(boxstyle='round', facecolor='white', alpha=0.5)

	    if flag_ecmwf:
		# place a text box in upper left in axes coords
		textstr = '$hBL=%.2f ~km$'%(out.hbl_ecmwf[0])
		if out.InJH_prmv2[0] > 0:
		    textstr+=' \n $\mathrm{PRM}v2=%.2f ~km$'%(out.InJH_prmv2[0])
		textstr2=''
		if out.InJH_prmv0[0] > 0:
		    textstr2+='$\mathrm{PRM}v0=%.2f ~km$'%(out.InJH_prmv0[0])
		if out.InJH_prmv1[0] > 0:
		    if textstr2 != '':
		        textstr2+='\n' 
		    textstr2+='$\mathrm{PRM}v1=%.2f ~km$'%(out.InJH_prmv1[0])
		if out.InJH_sof[0] > 0:
		    if textstr2 != '':
		        textstr2+='\n' 
		    textstr2+='$\mathrm{Sof}=%.2f ~km$'%(out.InJH_sof[0])

		ax.text(0.05, 0.12, textstr, transform=ax.transAxes, fontsize=24, verticalalignment='top', bbox=props)
		if textstr2 != '':
		    ax.text(0.3, 0.12, textstr2, transform=ax.transAxes, fontsize=24, verticalalignment='top', bbox=props)

	    print out.InJH_prmv0[0], out.InJH_prmv1[0], out.InJH_prmv2[0], out.InJH_sof[0]

	    fig.savefig(path_image+'geoRef_colorComposite_'+'{:04d}{:03d}{:04d}_{:s}.png'.format(year,doy,int(float(time)),sat))
	    plt.close(fig)
	    out.path_geoRef_ColorComp[0] = path_image+'geoRef_colorComposite_'+'{:04d}{:03d}{:04d}_{:s}.png'.format(year,doy,int(float(time)),sat)

    return out, grid_e, grid_n, dxy, conv_ll2utm, conv_utm2ll, AIcha, AIUncha, AISta, AIUnSta, Aero_geo, percent_obs
