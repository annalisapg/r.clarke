#!/usr/bin/env python
#%Module
#%  description: Calculates Clark hydrograph for a basin in output from r.traveltime.
#%  keywords: Clark, traveltime, runoff
#%End
#%option
#% key: dem
#% type: string
#% gisprompt: old,cell,raster
#% description: Dem of the zone
#% required: yes
#%end
#%option
#% key: manningsgrid
#% type: string
#% gisprompt: old,cell,raster
#% description: Basin roughness
#% required: yes
#%end
#%option
#% key: threshold
#% type: double
#% description: Threshold
#% required: yes
#%end
#%option
#% key: chanwidth
#% type: string
#% gisprompt: old,cell,raster
#% description: Channel width
#% required: yes
#%end
#%option
#% key: manningschan
#% type: string
#% gisprompt: old,cell,raster
#% description: Channel roughness
#% required: yes
#%end
#%option
#% key: adis
#% type: double
#% description: Average-year discharge in mc/s
#% required: yes
#%end
#%option
#% key: k
#% type: double
#% description: LAG of the basin (h)
#% required: yes
#%end
#%option
#% key: traveltime
#% type: string
#% gisprompt: new,cell,raster
#% description: Traveltime ouput map
#% required: yes
#%end
#%option
#% key: erain
#% type: string
#% gisprompt: old_file,file,input
#% key_desc: name
#% description: Text file of effective rain on the basin surface
#% required: yes
#%end
#%option
#% key: qtime
#% type: string
#% gisprompt: new_file,file,output
#% key_desc: name
#% description: Text file of q arriving to the cross section in time
#% required: yes
#%end
#%option
#% key: xout
#% type: double
#% description: xout
#% required: yes
#%end
#%option
#% key: yout
#% type: double
#% description: yout
#% required: yes
#%end

import sys
import os
import grass.script as grass
import math
from numpy import zeros
from numpy import array

def main():
	dem = options['dem'] 
	manningsgrid = options['manningsgrid']
	threshold = options['threshold']
	chanwidth = options['chanwidth']
	manningschan = options['manningschan']
	adis = options['adis']
	k = options['k']
	traveltime = options['traveltime']
	erain = options['erain']
	qtime = options['qtime']
	xout = options['xout']
	yout = options['yout']
#	grass.run_command('g.remove', rast = 'accu,rnetwork,drain,rnetwork_1,filled,filled_d' )
	grass.run_command('g.remove', rast = 'accu,rnetwork,drain,rnetwork_1,filled,filled_d,traveltime,traveltime_min' )
	grass.run_command('r.watershed', elevation = dem , threshold = threshold , accumulation = 'accu' , drainage = 'drain' , stream = 'rnetwork' )
	grass.run_command('r.null', map = 'rnetwork' , setnull = 0 )
	grass.run_command('r.mapcalculator' , amap = 'rnetwork' , formula = 'rnetwork/rnetwork' , outfile = 'rnetwork_1') 
	grass_region = grass.read_command('g.region' , flags = 'p').split()
#	g.region -p
	N = float(grass.read_command('g.region' , flags = 'ap').split('\n')[4].split(':')[-1])
	S = float(grass.read_command('g.region' , flags = 'ap').split('\n')[5].split(':')[-1])
	W = float(grass.read_command('g.region' , flags = 'ap').split('\n')[6].split(':')[-1])
	E = float(grass.read_command('g.region' , flags = 'ap').split('\n')[7].split(':')[-1])
	NSres_region = float(grass.read_command('g.region' , flags = 'ap').split('\n')[8].split(':')[-1])
	EWres_region = float(grass.read_command('g.region' , flags = 'ap').split('\n')[9].split(':')[-1])
	grass.run_command('g.region' , flags = 'ap', n = (N + NSres_region), s = (S - NSres_region) , w = (W - EWres_region) , e = (E + EWres_region) )
	grass.run_command('r.fill.dir', input = 'rnetwork', elevation = 'filled', direction = 'filled_d', type = 'grass')
	ladis = adis * 1000
	grass.run_command('r.traveltimeUp' , dir = 'drain', accu = 'accu', dtm = 'filled', manningsn = manningsgrid, out_x = xout, out_y = yout, threshold = threshold, nchannel = manningschan, b = chanwidth, dis = ladis, out = 'traveltime_min' )
	grass.run_command('r.mapcalculator' , amap = 'traveltime_min' , formula = 'traveltime_min/60' , outfile = 'traveltime')
	steps = int(float(grass.read_command('r.univar' , map = dem).split('\n')[7].split(':')[-1]))
	a = zeros((0,2), float)
	for i in grass.read_command('r.report', map = 'dem', units = 'k', null = '*', flags = 'h').replace('-','').replace('from  to . . . . . . . . . . . . . . . . .|','+').replace("|",'').split('\n'):
		try :
			i = float(i.split('+')[-1].strip())
			a = numpy.vstack([a, (i)])
		except :
			pass
	a = numpy.r_[[0], a]
	d = numpy.c_[numpy.arange(len(a)), a]
	e = zeros((len(d),2), float)
	f = zeros((len(d), 2), float)
	t = 1.026042
	lista = zeros((0, 4), float)
	cc = zeros((0,2), float)
	for i in range(len(d)):
		if d[i,0] == 0:
			e[i,0] , e[i,1] = d[i,0] , d[i,1]
			f[i,0] , f[i,1] = d[i,0] , d[i,1]
		if d[i,0] != 0:
			e[i,0] , e[i,1] = t * d[i,0] , d[i,1]
			f[i,0] , f[i,1] = (e[i,0] + e[i-1,0])/2 , d[i,1]
	for t in range(20*len(rain)+1):
		for i in range(len(f)):
			for j in range(len(rain)):
				if t >= (f[i,0] + rain[j,0]):
					fij=(f[i,1]*rain[j,1] * exp(- ( t - ( f[i,0] + rain[j,0] ) ) / 1.2 ) ) * 1000 / ( 1.2 * 3600 )
					lista = numpy.vstack([lista, (t, f[i,0], rain[j,0], fij)])
					
	for k in range(max(lista[:,0])+1):
		index = numpy.nonzero(lista[:,0] == k)[0]
		portata = sum(lista[index, 3])
		cc = numpy.vstack([cc, (k,portata)])
	plotImage(cc[:,0], cc[:,1],output_png,'-','Discharge [m^3/s]','Time [s]','Hydrograph (Clark)')

def plotImage(x,y,image,type,xlabel,ylabel,title):
	plt.plot(x, y, type)
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.title(title)
	plt.grid(True)
	plt.savefig(image)
	plt.close('all')


if __name__ == "__main__":
	options, flags = grass.parser()
	sys.exit(main())

