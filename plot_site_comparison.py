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
ano=2008
sat=1


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
satelites='VGT', 'ATS', 'AT2'
site_alt=['CAN', 'COL', 'BRA', 'POR', 'ANG', 'SAF', 'KAZ', 'BOR', 'RUS', 'AUS']
year=str(ano)
year2=year

sizex=18
sizey=16



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



meses='01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'
#File = '/home/bmota/CCI_fire/Results/Version_12/VGT/VGT_por_SS04_2005_v12_Nov12/BA_PIX_VGT_SS' #'01.tif'

File = '/home/bmota/CCI_fire/Results/Version_12/VGT/VGT_'+sitest[site-1]+'_SS'+stite[site-1]+'_'+year+'_v12_Nov12/BA_PIX_VGT_SS' #'01.tif'
File2 = '/home/bmota/CCI_fire/Results/Global/Site_Comparison/CopySS_2008_1km'+site_alt[site-1]+'.tif'
# le os dados
print(str(File+stite[site-1]+'_'+year+'.tif'))
dataset = gdal.Open(str(File+stite[site-1]+'_'+year2+'.tif'), GA_ReadOnly )
data_1 = dataset.GetRasterBand(1).ReadAsArray()
data_2 = dataset.GetRasterBand(2).ReadAsArray()
data_3 = dataset.GetRasterBand(3).ReadAsArray()
data_4 = dataset.GetRasterBand(4).ReadAsArray()
data_5 = dataset.GetRasterBand(5).ReadAsArray()
data_6 = dataset.GetRasterBand(6).ReadAsArray()

data_1a=data_1
data_2a=data_2
data_3a=data_3
data_4a=data_4
data_5a=data_5
data_6a=data_6

#mascarra os dados para valores invalidos
data_1 = np.ma.masked_equal(data_1a, 0)
data_1 = np.ma.masked_equal(data_1, -12)
data_2 = np.ma.masked_equal(data_2a, 0)
data_2 = np.ma.masked_equal(data_2, -12)
data_3 = np.ma.masked_equal(data_3a, 0)
data_3 = np.ma.masked_equal(data_3, -12)
data_4 = np.ma.masked_equal(data_4a, 0)
data_4 = np.ma.masked_equal(data_4, -12)
data_5 = np.ma.masked_equal(data_5a, 0)
data_5 = np.ma.masked_equal(data_5, -12)
data_6 = np.ma.masked_equal(data_6a, 0)
data_6 = np.ma.masked_equal(data_6, -12)

dataset2 = gdal.Open(str(File2), GA_ReadOnly )
data_1b = dataset2.GetRasterBand(1).ReadAsArray()
data_2b = dataset2.GetRasterBand(2).ReadAsArray()
data_3b = dataset2.GetRasterBand(3).ReadAsArray()
data_4b = dataset2.GetRasterBand(4).ReadAsArray()
data_5b = dataset2.GetRasterBand(5).ReadAsArray()
data_6b = dataset2.GetRasterBand(6).ReadAsArray()

data_1c=data_1b
data_2c=data_2b
data_3c=data_3b
data_4c=data_4b
data_5c=data_5b
data_6c=data_6b

#mascarra os dados para valores invalidos
data_1b = np.ma.masked_equal(data_1c, 0)
data_1b = np.ma.masked_equal(data_1b, -12)
data_2b = np.ma.masked_equal(data_2c, 0)
data_2b = np.ma.masked_equal(data_2b, -12)
data_3b = np.ma.masked_equal(data_3c, 0)
data_3b = np.ma.masked_equal(data_3b, -12)
data_4b = np.ma.masked_equal(data_4c, 0)
data_4b = np.ma.masked_equal(data_4b, -12)
data_5b = np.ma.masked_equal(data_5c, 0)
data_5b = np.ma.masked_equal(data_5b, -12)
data_6b = np.ma.masked_equal(data_6c, 0)
data_6b = np.ma.masked_equal(data_6b, -12)







# ------------------------------------------------------------------------------------------------------
# Banda 1 - dia de deteccao
# ------------------------------------------------------------------------------------------------------
if tag==1:
	fig=plt.figure()
else:
	fig=plt.figure(num=None, figsize=(sizex, sizey), dpi=300, facecolor='w', edgecolor='k')

band=1
data=data_1

fig.add_subplot(121)
#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=366, cmap=plt.cm.hsv)
titulo='day'
titulo_final='Day of burnt'+' of '+siteslong[site-1]+' site for year '+year

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


data=data_1b

fig.add_subplot(122)
#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=366, cmap=plt.cm.hsv)
titulo='day'
titulo_final='Day of burnt'+' of '+siteslong[site-1]+' site for year '+year

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

if tag==1:
	plt.show()
else:
	plt.savefig(str(satelites[sat-1]+'_Comp_SS'+stite[site-1]+'_'+str(year)+'_band'+str(band)+'.png'), dpi=300)


# ------------------------------------------------------------------------------------------------------
# Banda 2 - confidence level
# ------------------------------------------------------------------------------------------------------
if tag==1:
	fig=plt.figure()
else:
	fig=plt.figure(num=None, figsize=(sizex, sizey), dpi=300, facecolor='w', edgecolor='k')


fig.add_subplot(121)
band=2
data=data_2

#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=100, cmap=plt.cm.jet)
titulo='conf.'
titulo_final='confidence interval'+' of '+siteslong[site-1]+' site for year '+year

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

data=data_2b

fig.add_subplot(122)
#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=100, cmap=plt.cm.jet)
titulo='conf.'
titulo_final='confidence interval'+' of '+siteslong[site-1]+' site for year '+year

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

if tag==1:
	plt.show()
else:
	plt.savefig(str(satelites[sat-1]+'_Comp_SS'+stite[site-1]+'_'+str(year)+'_band'+str(band)+'.png'), dpi=300)


# ------------------------------------------------------------------------------------------------------
# Banda 3 - Lag day
# ------------------------------------------------------------------------------------------------------
if tag==1:
	fig=plt.figure()
else:
	fig=plt.figure(num=None, figsize=(sizex, sizey), dpi=300, facecolor='w', edgecolor='k')



fig.add_subplot(121)
band=3
data=data_3

#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=15, cmap=plt.cm.jet)
titulo='Days'
titulo_final='Day lag of detection'+' of '+siteslong[site-1]+' site for year '+year

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

data=data_3b

fig.add_subplot(122)
#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=15, cmap=plt.cm.jet)
titulo='Days'
titulo_final='Day lag of detection'+' of '+siteslong[site-1]+' site for year '+year

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

if tag==1:
	plt.show()
else:
	plt.savefig(str(satelites[sat-1]+'_Comp_SS'+stite[site-1]+'_'+str(year)+'_band'+str(band)+'.png'), dpi=300)


# ------------------------------------------------------------------------------------------------------
# Banda 4 - Pixel day coverage
# ------------------------------------------------------------------------------------------------------
if tag==1:
	fig=plt.figure()
else:
	fig=plt.figure(num=None, figsize=(sizex, sizey), dpi=300, facecolor='w', edgecolor='k')


fig.add_subplot(121)
band=4
data=data_4

#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
if sat==1:
	im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=366, cmap=plt.cm.jet)
else:
	im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=60, cmap=plt.cm.jet)

titulo='obs.'
titulo_final='Valid pixels'+' of '+siteslong[site-1]+' site for year '+year



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
data=data_4b

fig.add_subplot(122)
#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=366, cmap=plt.cm.jet)
titulo='obs'
titulo_final='Valid pixels'+' of '+siteslong[site-1]+' site for year '+year

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

if tag==1:
	plt.show()
else:
	plt.savefig(str(satelites[sat-1]+'_Comp_SS'+stite[site-1]+'_'+str(year)+'_band'+str(band)+'.png'), dpi=300)


# ------------------------------------------------------------------------------------------------------
# Banda 5 - Valid surface observations
# ------------------------------------------------------------------------------------------------------
if tag==1:
	fig=plt.figure()
else:
	fig=plt.figure(num=None, figsize=(sizex, sizey), dpi=300, facecolor='w', edgecolor='k')


fig.add_subplot(121)
band=5
data=data_5

#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
if sat==1:
	im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=366, cmap=plt.cm.jet)
else:
	im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=150, cmap=plt.cm.jet)
titulo='obs.'
titulo_final='Pixel coverage'+' of '+siteslong[site-1]+' site for year '+year



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


data=data_5b

fig.add_subplot(122)
#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=366, cmap=plt.cm.jet)
titulo='obs'
titulo_final='Pixel coverage'+' of '+siteslong[site-1]+' site for year '+year

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

if tag==1:
	plt.show()
else:
	plt.savefig(str(satelites[sat-1]+'_Comp_SS'+stite[site-1]+'_'+str(year)+'_band'+str(band)+'.png'), dpi=300)


# ------------------------------------------------------------------------------------------------------
# Banda 6 - Cloudy observations
# ------------------------------------------------------------------------------------------------------
if tag==1:
	fig=plt.figure()
else:
	fig=plt.figure(num=None, figsize=(sizex, sizey), dpi=300, facecolor='w', edgecolor='k')


fig.add_subplot(121)
band=6
data=data_6

#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=100, cmap=plt.cm.jet)
titulo='obs.'
titulo_final='Cloud pixels'+' of '+siteslong[site-1]+' site for year '+year


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
#if tag==1:
#	plt.show()
#else:
#	plt.savefig(str(siteslong[site-1]+year+'map.png'), dpi=600)
#savefig('map.png', dpi=600)

data=data_6b

fig.add_subplot(122)
#estabelece a projeccao
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[1],urcrnrlat=lat[1],resolution='i')

#obtencao das coordenadas dos cantos e definicao dos extents
canto=dataset.GetGeoTransform()
nLat, nLon  = data.shape
LL = p1(canto[0], canto[3]+canto[5]*data.shape[0], inverse=True)
UR = p1(canto[0]+canto[1]*data.shape[1], canto[3], inverse=True)

# desenha mascarra da agua
m.drawlsmask(land_color='0.3', ocean_color='w', lakes=True, resolution='f')

#calculo dos extents do mapa projectado
extent = (LL[0], UR[0], LL[1], UR[1])
x, y = m(extent[:2],extent[2:])
extent_xy = (x[0], x[1], y[0], y[1])

# plot da imagem
im=plt.imshow(data, extent=extent_xy, vmin=1, vmax=100, cmap=plt.cm.jet)
titulo='obs.'
titulo_final='Cloud pixels'+' of '+siteslong[site-1]+' site for year '+year

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

if tag==1:
	plt.show()
else:
	plt.savefig(str(satelites[sat-1]+'_Comp_SS'+stite[site-1]+'_'+str(year)+'_band'+str(band)+'.png'), dpi=300)



