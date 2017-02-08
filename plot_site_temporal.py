from __future__ import print_function
from mpl_toolkits.basemap import pyproj
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon 
import osgeo.gdal as gdal
from osgeo.gdalconst import *
import numpy as np
import matplotlib as mpl
import sys


#sat=1
#band=1
#version=12
tag=np.int(sys.argv[1])
site=np.int(sys.argv[2])
#band=1
band=np.int(sys.argv[3])
version=np.int(sys.argv[4])
sat=np.int(sys.argv[5])

sizex=18
sizey=16

if tag==1:
	mpl.rcParams['font.size'] = 10.
	#mpl.rcParams['font.family'] = 'Comic Sans MS'
	mpl.rcParams['axes.labelsize'] = 8.
	mpl.rcParams['xtick.labelsize'] = 6.
	mpl.rcParams['ytick.labelsize'] = 6.
	font=8
else:
	mpl.rcParams['font.size'] = 6.
	#mpl.rcParams['font.family'] = 'Comic Sans MS'
	mpl.rcParams['axes.labelsize'] = 4.
	mpl.rcParams['xtick.labelsize'] = 4.
	mpl.rcParams['ytick.labelsize'] = 4.
	font=4

stite=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10']
sitest=['can', 'col', 'bra', 'por', 'ang', 'saf', 'kaz', 'bor', 'rus', 'aus']
siteslong=['Canada', 'Colombia', 'Brazil', 'Portugal', 'Angola', 'South Africa', 'Kazakistan', 'Borneo', 'Russia', 'Australia']
meses='01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'
satelites='VGT', 'ATS', 'AT2'

if site==1:
	p1 = pyproj.Proj('+proj=utm +zone=13 +datum=WGS84')
	lon_0=-105.
	lat_0=0.
	lat=52., 63.
	lon=-112., -104.
if site==2:
	p1 = pyproj.Proj('+proj=utm +zone=19 +datum=WGS84')
	lon_0=-69.
	lat_0=0.
	lat=3., 9.
	lon=-73., -66.
if site==3:
	p1 = pyproj.Proj('+proj=utm +zone=22 +south +datum=WGS84')
	lon_0=-51.
	lat_0=0.
	lat=-13., -6.
	lon=-52., -46.
if site==4:
	p1 = pyproj.Proj('+proj=utm +zone=29 +datum=WGS84')
	lon_0=-9.
	lat_0=0.
	lat=37., 45.
	lon=-10., -4.
if site==5:
	p1 = pyproj.Proj('+proj=utm +zone=33 +south +datum=WGS84')
	lon_0=15.
	lat_0=0.
	lat=-12., -5.
	lon=13, 19.
if site==6:
	p1 = pyproj.Proj('+proj=utm +zone=35 +south +datum=WGS84')
	lon_0=27.
	lat_0=0.
	lat=-28., -21.
	lon=27, 33.
if site==7:
	p1 = pyproj.Proj('+proj=utm +zone=40 +datum=WGS84')
	lon_0=57.
	lat_0=0.
	lat=46., 54.
	lon=53., 59.
if site==8:
	p1 = pyproj.Proj('+proj=utm +zone=49 +south +datum=WGS84')
	lon_0=111.
	lat_0=0.
	lat=-5., 1.
	lon=111, 117.
if site==9:
	p1 = pyproj.Proj('+proj=utm +zone=49 +datum=WGS84')
	lon_0=111.
	lat_0=0.
	lat=49., 59.
	lon=109, 117.
if site==10:
	p1 = pyproj.Proj('+proj=utm +zone=53 +south +datum=WGS84')
	lon_0=135.
	lat_0=0.
	lat=-18., -11.
	lon=129., 136.




if sat==1:
	years1=range(1998, 2004);
	years2=range(2004, 2010);

if sat==2: 
	years1=range(2002, 2008);
	years2=range(2008, 2009);
if sat==3: 
	years1=range(1995, 2001);
	years2=range(2001, 2003);

fig=plt.figure(num=None, figsize=(sizex, sizey), dpi=600, facecolor='w', edgecolor='k')

ind=0
for ano in years1:
	ind=ind+1
	print(ind)	

	fig.add_subplot(2,3,ind)

	year=str(ano)
#	if sat==1 and version==12:
#		if site==4 or site==5 or site==6 or site==8 or site==9:
#			if ano<2000:
#				year2=str(ano-1900)
#			else:
#				year2=str(ano-2000)
#		else:
#			year2=year
#	else:
	year2=year
	if version==11:
		File = '/home/bmota/CCI_fire/Results/Version_11/'+satelites[sat-1]+'/'+satelites[sat-1]+'_'+sitest[site-1]+'_SS'+stite[site-1]+'_'+year+'_v11_Aug12/BA_PIX_'+satelites[sat-1]+'_SS' 
	if version==12:
		File = '/home/bmota/CCI_fire/Results/Version_12/'+satelites[sat-1]+'/'+satelites[sat-1]+'_'+sitest[site-1]+'_SS'+stite[site-1]+'_'+year+'_v12_Nov12/BA_PIX_'+satelites[sat-1]+'_SS' 

	print(str(File+stite[site-1]+'_'+year+meses[0]+'.tif'))
	dataset = gdal.Open(str(File+stite[site-1]+'_'+year2+meses[0]+'.tif'), GA_ReadOnly )
	data_1 = dataset.GetRasterBand(band).ReadAsArray()
	data_1a=data_1*1.0

	mes=range(2,13)
	for i in mes:
		print(str(File+stite[site-1]+'_'+year+meses[i-1]+'.tif'))
		dataset = gdal.Open(str(File+stite[site-1]+'_'+year2+meses[i-1]+'.tif'), GA_ReadOnly )
		data_1 = dataset.GetRasterBand(band).ReadAsArray()
		data_1a = data_1a + data_1

	#mascarra os dados para valores invalidos
	data_1 = np.ma.masked_equal(data_1a, 0)
	data_1 = np.ma.masked_equal(data_1, -12)

	data=data_1

	#estabelece a projeccao
	m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

	#obtencao das coordenadas dos cantos e definicao dos extents
	canto=dataset.GetGeoTransform()
	#nLat, nLon  = data.shape
	LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
	UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

	# desenha mascarra da agua
	m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

	#calculo dos extents do mapa projectado
	extent = (LL[0], UR[0], LL[1], UR[1])
	x, y = m(extent[:2],extent[2:])
	extent_xy = (x[0], x[1], y[0], y[1])

	# plot da imagem
	if band==1:
		im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=366, cmap=plt.cm.hsv)
		titulo='day'
		titulo_final='Day of burnt'+' of '+siteslong[site-1]+' site for year '+year
	if band==2:
		im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=100, cmap=plt.cm.jet)
		titulo='conf.layer'
		titulo_final='confidence Interval'+' of '+siteslong[site-1]+' site for year '+year
	if band==3:
		im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=15, cmap=plt.cm.jet)
		titulo='Day lag'
		titulo_final='Day lag of detection'+' of '+siteslong[site-1]+' site for year '+year

	if band==4 or band==5:
		im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=366, cmap=plt.cm.jet)
		titulo='obs.'
		if band==4:
			titulo_final='Valid observation'+' of '+siteslong[site-1]+' site for year '+year
		if band==5:
			titulo_final='Covered observations'+' of '+siteslong[site-1]+' site for year '+year
	if band==6:
		im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=100, cmap=plt.cm.jet)
		titulo='obs.'
		titulo_final='Cloud observations'+' of '+siteslong[site-1]+' site for year '+year


	#plot do frame a volta do study site
	p = Polygon([(x[0], y[0]), (x[1], y[0]), (x[1], y[1]), (x[0], y[1])],facecolor='none',edgecolor='black',linewidth=2)
	plt.gca().add_patch(p) 

	#plot da barra de
	cbar = m.colorbar(im,"right", size="5%",pad="2%")
	c1=cbar.set_label(titulo)
	for t in cbar.ax.get_yticklabels():
     		t.set_fontsize(font)

	#caracteristicas do mapa
	ax=m.ax
	m.drawparallels(range(int(np.array(lat).min()),int(np.array(lat).max())),labels=[1,0,0,0],labelstyle='+/-')
	m.drawmeridians(range(int(np.array(lon).min()),int(np.array(lon).max())),labels=[0,0,0,1],labelstyle='+/-')
	m.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
	m.drawcoastlines(linewidth=0.5)
	m.drawcountries(linewidth=0.5)

	#titulo do mapa
	plt.title(titulo_final)



#plt.ylabel("Latitude")
#plt.xlabel("Longitude")
if tag==1:
	plt.show()
else:
	plt.savefig(str(satelites[sat-1])+'_'+str(siteslong[site-1]+'_band_'+str(band)+'_v'+str(version)+'_part1map.png'), dpi=600)
#savefig('map.png', dpi=600)

fig=plt.figure(num=None, figsize=(sizex, sizey), dpi=600, facecolor='w', edgecolor='k')

ind=0
for ano in years2:
	ind=ind+1
	print(ind)	

	fig.add_subplot(2,3,ind)

	year=str(ano)
#	if sat==1 and version==12:
#		if site==4 or site==5 or site==6 or site==8 or site==9:
#			if ano<2000:
#				year2=str(ano-1900)
#			else:
#				year2=str(ano-2000)
#		else:
#			year2=year
#	else:
	year2=year

	if version==11:
		File = '/home/bmota/CCI_fire/Results/Version_11/'+satelites[sat-1]+'/'+satelites[sat-1]+'_'+sitest[site-1]+'_SS'+stite[site-1]+'_'+year+'_v11_Aug12/BA_PIX_'+satelites[sat-1]+'_SS' 
	if version==12:
		File = '/home/bmota/CCI_fire/Results/Version_12/'+satelites[sat-1]+'/'+satelites[sat-1]+'_'+sitest[site-1]+'_SS'+stite[site-1]+'_'+year+'_v12_Nov12/BA_PIX_'+satelites[sat-1]+'_SS' 

	print(str(File+stite[site-1]+'_'+year+meses[0]+'.tif'))
	dataset = gdal.Open(str(File+stite[site-1]+'_'+year2+meses[0]+'.tif'), GA_ReadOnly )
	data_1 = dataset.GetRasterBand(band).ReadAsArray()
	data_1a=data_1*1.0

	mes=range(2,13)
	for i in mes:
		print(str(File+stite[site-1]+'_'+year+meses[i-1]+'.tif'))
		dataset = gdal.Open(str(File+stite[site-1]+'_'+year2+meses[i-1]+'.tif'), GA_ReadOnly )
		data_1 = dataset.GetRasterBand(band).ReadAsArray()
		data_1a = data_1a + data_1

	#mascarra os dados para valores invalidos
	data_1 = np.ma.masked_equal(data_1a, 0)
	data_1 = np.ma.masked_equal(data_1, -12)

	data=data_1

	#estabelece a projeccao
	m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

	#obtencao das coordenadas dos cantos e definicao dos extents
	canto=dataset.GetGeoTransform()
	#nLat, nLon  = data.shape
	LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
	UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

	# desenha mascarra da agua
	m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

	#calculo dos extents do mapa projectado
	extent = (LL[0], UR[0], LL[1], UR[1])
	x, y = m(extent[:2],extent[2:])
	extent_xy = (x[0], x[1], y[0], y[1])

	# plot da imagem
	if band==1:
		im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=366, cmap=plt.cm.hsv)
		titulo='day'
		titulo_final='Day of burnt'+' of '+siteslong[site-1]+' site for year '+year
	if band==2:
		im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=100, cmap=plt.cm.jet)
		titulo='conf.layer'
		titulo_final='confidence Interval'+' of '+siteslong[site-1]+' site for year '+year
	if band==3:
		im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=15, cmap=plt.cm.jet)
		titulo='Day lag'
		titulo_final='Day lag of detection'+' of '+siteslong[site-1]+' site for year '+year

	if band==4 or band==5:
		im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=366, cmap=plt.cm.jet)
		titulo='obs.'
		if band==4:
			titulo_final='Valid observation'+' of '+siteslong[site-1]+' site for year '+year
		if band==5:
			titulo_final='Covered observations'+' of '+siteslong[site-1]+' site for year '+year
	if band==6:
		im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=100, cmap=plt.cm.jet)
		titulo='obs.'
		titulo_final='Cloud observations'+' of '+siteslong[site-1]+' site for year '+year


	#plot do frame a volta do study site
	p = Polygon([(x[0], y[0]), (x[1], y[0]), (x[1], y[1]), (x[0], y[1])],facecolor='none',edgecolor='black',linewidth=2)
	plt.gca().add_patch(p) 

	#plot da barra de
	cbar = m.colorbar(im,"right", size="5%",pad="2%")
	c1=cbar.set_label(titulo)
	for t in cbar.ax.get_yticklabels():
     		t.set_fontsize(font)

	#caracteristicas do mapa
	ax=m.ax
	m.drawparallels(range(int(np.array(lat).min()),int(np.array(lat).max())),labels=[1,0,0,0],labelstyle='+/-')
	m.drawmeridians(range(int(np.array(lon).min()),int(np.array(lon).max())),labels=[0,0,0,1],labelstyle='+/-')
	m.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
	m.drawcoastlines(linewidth=0.5)
	m.drawcountries(linewidth=0.5)

	#titulo do mapa
	plt.title(titulo_final)



#plt.ylabel("Latitude")
#plt.xlabel("Longitude")
if tag==1:
	plt.show()
else:
	plt.savefig(str(satelites[sat-1])+'_'+str(siteslong[site-1]+'_band_'+str(band)+'_v'+str(version)+'_part2map.png'), dpi=600)
#savefig('map.png', dpi=600)
