#!/usr/bin/env python
# coding: utf-8

#modules utiles pour python
# matplotlib = pour faire les figures
# numpy : trucs numériques
# pandas: dataframes gestion "excel" 
# os: pour changer de directory, ...
import matplotlib.pyplot as plt  #plt.FUNCTIONNAE
import numpy as np  # outils numériques
import pandas as pd  # gestion csv, dataframe
import os # pour changer de directory


def main():
	is1D_x = False 	#1s1D_x = False # si c'est 1D en y
	plotAcidity = False
	generic_title = "Enviro_" # pour mettre un titre global
	profile_sample = [118]  # si tu veux aucun profile_sample, faut écrire: profile_sample = []
	folder = []
	#folder.append(Simulation("TSU13_Prodata", legend_label="Enviro_Prodata",color="orange"))
	folder.append(Simulation("TSU13_Chess", legend_label="Enviro_Chess",color="b"))

	node_list = []  # [] si tu veux pas de node_list # numero du noeud
	fnode_list = []  # [] si tu veux pas de node_list # numero du noeud
	lag_time = 0 # tjs en année
	ListOfVariables = []

	if(is1D_x):
		ListOfVariables.append(plot_var("x-distance",unit = "(m)",fancy_name=" x (m)",to_plot=False,transient=False))
	else:
		ListOfVariables.append(plot_var("y-distance",unit = "(m)",fancy_name=" y (m)",to_plot=False,transient=False))

	ListOfVariables.append(plot_var("pH", ymin=1, ymax=9))  # tmin et tmax = sample
#	ListOfVariables.append(plot_var("Eh"))
#	ListOfVariables.append(plot_var("aqueous{Ca[2+]}",fancy_name="Ca(aq)", fname = "Ca(aq)",type_sort=1,color="r"))
#	ListOfVariables.append(plot_var("aqueous{Mg[2+]}",fancy_name="Mg(aq)", fname = "Mg(aq)",type_sort=1,color="g"))
	ListOfVariables.append(plot_var("aqueous{UO2[2+]}",fancy_name="UO2(aq)", fname = "UO2(aq)", logscale= True, type_sort=1,color="b"))
	ListOfVariables.append(plot_var("aqueous{SO4[2-]}",fancy_name="SO4(aq)", fname = "SO4(aq)", logscale= True, type_sort=1,color="orange"))
#	ListOfVariables.append(plot_var("aqueous{K[+]}",fancy_name="K(aq)", fname = "K(aq)", type_sort=1,color="c"))
#	ListOfVariables.append(plot_var("aqueous{Na[+]}", fancy_name="Na(aq)", fname = "Na(aq)", type_sort=1,color="m"))
	ListOfVariables.append(plot_var("fixed{H[+]}", fancy_name="Fixed H+", fname = "Fixed H+",  type_sort=2,color="k"))
	ListOfVariables.append(plot_var("fixed{Ca[2+]}",fancy_name="Fixed Ca++", fname = "Fixed Ca++",type_sort=2,color="r"))
	ListOfVariables.append(plot_var("fixed{Mg[2+]}",fancy_name="Fixed Mg++", fname = "Fixed Mg++",type_sort=2,color="g"))
	ListOfVariables.append(plot_var("fixed{UO2[2+]}",fancy_name="Fixed UO2++", fname = "Fixed UO2++", type_sort=2,color="b"))
	ListOfVariables.append(plot_var("fixed{K[+]}",fancy_name="Fixed K+", fname = "Fixed K+",type_sort=2,color="c"))
	ListOfVariables.append(plot_var("fixed{Na[+]}",fancy_name="Fixed Na+", fname = "Fixed Na+",type_sort=2,color="m"))
	

	Combined_plot_properties = []
	Combined_plot_properties.append(Combined_plot("Aqueous","aq","Aqueous (mol/L)",logscale=True))
	Combined_plot_properties.append(Combined_plot("Sorbed","fix","Sorbed (mol/L)",logscale=True))

	FluxVariables = []
	FluxVariables.append(plot_var("aqueous{UO2[2+]}", fname = "U(t)", unit = "T", fac = 0.000238,tmin=0, tmax=50, ymin=0, ymax=2))
	if(plotAcidity):
		FluxVariables.append(plot_var("acidity",acidity=True,unit = "kg",tmin=0,tmax=50,ymin=0,ymax=2))
        
	transient_data = []
	profile_data = []
	prop_cycle = plt.rcParams['axes.prop_cycle']
	colors = prop_cycle.by_key()['color']
	#print("number of colors in cycle is %d" %(len(colors)))
	for i in range(len(ListOfVariables)):
	    if(ListOfVariables[i].color==""):
	        ListOfVariables[i].color = "b"
	    if(ListOfVariables[i].transient):
	        transient_data.append(ListOfVariables[i])            
	    if(ListOfVariables[i].profile):
	        profile_data.append(ListOfVariables[i])

	for i in range(len(folder)):
	    if(folder[i].color==""):
	        folder[i].color=colors[i]
	for i in range(len(FluxVariables)):
	    if(FluxVariables[i].color==""):
	        FluxVariables[i].color=colors[i]
	Profile_res = []
	Transient_res = []
	profile_data_names = []
	transient_data_names = []

	for i in range(len(profile_data)):
		profile_data_names.append(profile_data[i].name)
    
	for i in range(len(transient_data)):
		transient_data_names.append(transient_data[i].name)
    
	NumberOfSampleTimes = 0
	for i in range(len(folder)):
		fold = folder[i].fold
		varNames, sTimes, Nnodes, FirstLine, units = extractData(fold)
		if(i==0):
			Sample_times = np.zeros((len(folder),len(sTimes)))
            
		if(len(sTimes)>NumberOfSampleTimes):
 			if(i==0):
				 Sample_times = np.zeros((len(folder),len(sTimes)))
 			else:
				 Sample_times = extend_vector(Sample_times,len(sTimes),len(folder),NumberOfSampleTimes)
 			NumberOfSampleTimes = len(sTimes)
    
		#print(sTimes)
		#print(len(Sample_times))
		#print(len(Sample_times[0]))
		for k in range(len(sTimes)):
 			Sample_times[i][k] = sTimes[k]

		for sam in range(len(profile_sample)):
			Res = get_profile_data(profile_data_names,profile_sample[sam],fold,FirstLine,Nnodes,varNames)
			Profile_res.append(Res)
        
    	# 2. get data for transient plots
		for node in range(len(node_list)):
			Res = get_transient_data(transient_data_names,fold,FirstLine,node_list[node],Nnodes,varNames,sTimes)  
			Transient_res.append(Res)


	render_time = []
	t_unit = []
	for i in range(len(sTimes)):
		if(sTimes[i]<3600*24):
		    render_time.append(int(sTimes[i])-lag_time*365*3600*24)
		    if(render_time[len(render_time)-1]==1):
		        t_unit.append("second")
		    else:
		        t_unit.append("seconds")
		elif(sTimes[i]<3600*24*30):
		    render_time.append(int(sTimes[i]/3600/24)-lag_time*365)
		    if(render_time[len(render_time)-1]==1):
		        t_unit.append("day")
		    else:
		        t_unit.append("days")
		elif(sTimes[i]<3600*24*365):
		    render_time.append(int(sTimes[i]/3600/24/30)-lag_time*12)
		    if(render_time[len(render_time)-1]==1):
		        t_unit.append("month")
		    else:          
		        t_unit.append("months")
		else:
		    render_time.append(int(sTimes[i]/3600/24/365)-lag_time)
		    if(render_time[len(render_time)-1]==1):
		        t_unit.append("year")
		    else:
		        t_unit.append("years")

	plot_t_unit = []
	#plot_time = []
	if(sTimes[len(sTimes)-1]<3600*24):
		fac = 1
		plot_t_unit = "seconds"
	elif(sTimes[len(sTimes)-1]<3600*24):
		fac = 3600*24
		plot_t_unit = "days"
	elif(sTimes[len(sTimes)-1]<3600*24*365):
		fac = 3600*24*12
		plot_t_unit = "months"
	elif(sTimes[len(sTimes)-1]>=3600*24*365):
		fac = 3600*24*365
		plot_t_unit = "years"
		
	for j in range(len(folder)):
		sTimes= Sample_times[j][:]
		for i in range(len(sTimes)):
			if(sTimes[i]==0 and i>1):
				i=len(sTimes)                
			else:
				folder[j].plot_time.append((sTimes[i]-lag_time/3600/24/365)/fac)
		

	for i in range(len(profile_sample)):
		for j in range(len(profile_data)):
		    if(profile_data[j].to_plot):
		        fig, ax = plt.subplots()
		        for f in range(len(folder)):
		            df = Profile_res[i+f*(len(profile_sample))] 
		            if(is1D_x):
		                xval = df["x-distance"]
		                xlab = "x-distance (m)"
		            else:
		                xval = df["y-distance"]
		                xlab = "y-distance (m)"
		            if(len(folder)>1):
		                plot_color = folder[f].color
		            else:
		                plot_color = profile_data[j].color
		            if(profile_data[j].logscale):
		                if(profile_data[j].name in df.columns):
		                    plt.semilogy(xval,profile_data[j].fac*df[profile_data[j].name], label=str(folder[f].legend_label),color=plot_color)
		            else:
		                if(profile_data[j].name in df.columns):
		                    plt.plot(xval,profile_data[j].fac*df[profile_data[j].name], label=str(folder[f].legend_label),color=plot_color)
                            
		        if(profile_data[j].setxmax==True):
		            ax.set_xlim(right=profile_data[j].xmax)
		        if(profile_data[j].setxmin==True):
		            ax.set_xlim(left=profile_data[j].xmin)
		        if(profile_data[j].setymax==True):
		            ax.set_ylim(top=profile_data[j].ymax)
		        if(profile_data[j].setxmin==True):
		            ax.set_ylim(bottom=profile_data[j].ymin)

		        title = str(profile_data[j].fancy_name)+" after "+str(render_time[profile_sample[i]])+" "+t_unit[profile_sample[i]]        
		        plt.xlabel(xlab)
		        ylab = profile_data[j].fancy_name+" ("+profile_data[j].unit+")" #= "Time ("+plot_t_unit+")"
		        plt.ylabel(ylab)
		        plt.suptitle(title)
		        ax.legend()

		        filename = generic_title+str(profile_data[j].fname)+"_sample"+str(profile_sample[i])+".png"
		        fname = ''.join(filename)
		        #print(fname)
		        if(len(folder)==1):
		            #os.chdir(folder[f].fold)
		            plt.savefig(str(filename),dpi=500) #,transparent=True)
		            #os.chdir("..")
		        else:
		            plt.savefig(str(filename),dpi=500) #,transparent=True)
		        plt.close()
         

#################### TRANSIENT PLOTS #################################
	for i in range(len(node_list)):
		for j in range(len(transient_data)):
		    fig, ax = plt.subplots()
		    for f in range(len(folder)):
		        df = Transient_res[i+f*len(node_list)]  
		        if(len(folder)>1):
		            plot_color = folder[f].color
		        else:
		            plot_color = transient_data[j].color
		        if(transient_data[j].logscale):
		            if(transient_data[j].name in df.columns):
		                plt.semilogy(folder[f].plot_time,transient_data[j].fac*df[transient_data[j].name], label=str(folder[f].legend_label),color=plot_color)
		        else:
		            if(transient_data[j].name in df.columns):
		                plt.plot(folder[f].plot_time,transient_data[j].fac*df[transient_data[j].name], label=str(folder[f].legend_label),color=plot_color)

		    if(transient_data[j].settmax==True):
		        ax.set_xlim(right=folder[f].plot_time[transient_data[j].tmax])
		    if(transient_data[j].settmin==True):
		        ax.set_xlim(left=folder[f].plot_time[transient_data[j].tmin])

		    if(transient_data[j].setymax==True):
		        ax.set_ylim(top=transient_data[j].ymax)
		    if(transient_data[j].setymin==True):
		        ax.set_ylim(bottom=transient_data[j].ymin)


		    title = str("Transient "+transient_data[j].fancy_name+" at node "+str(node_list[i]))   
		    xlab = "Time ("+plot_t_unit+")"
		    plt.xlabel(xlab)
		    ylab = transient_data[j].fancy_name+" ("+transient_data[j].unit+")" 
		    plt.ylabel(ylab)
		    plt.suptitle(title)
		    ax.legend()

		    filename = generic_title+str(transient_data[j].fname)+"_node"+str(node_list[i])+".png"
		    fname = ''.join(filename)
		    #print(fname)
		    if(len(folder)==1):
		        #os.chdir(folder[f].fold)
		        plt.savefig(str(filename),dpi=500) #,transparent=True)
		        #plt.clf()
		        #os.chdir("..")
		    else:          
		        plt.savefig(str(filename),dpi=500) #,transparent=True)
		plt.close()


	######################## COMBINED PLOTS ############################
	max_sort = 0
	for i in range(len(profile_data)):
		if(profile_data[i].type_sort>max_sort):
		    max_sort = profile_data[i].type_sort

	for f in range(len(folder)):
		for sort in range(1,max_sort+1):
		    for i in range(len(profile_sample)):
		        #for j in range(1,max_sort): 
		        fig, ax = plt.subplots()        
		        df = Profile_res[i+f*(len(profile_sample))] 
		        if(is1D_x):
		            xval = df["x-distance"]
		            xlab = "x-distance (m)"
		        else:
		            xval = df["y-distance"]
		            xlab = "y-distance (m)"
		        
		        for j in range(len(profile_data)):
		            if(profile_data[j].type_sort==sort):                            
		                plot_color=profile_data[j].color
		                if(profile_data[j].name in df):
		                    if(Combined_plot_properties[sort-1].logscale == True):
		                        plt.semilogy(xval,profile_data[j].fac*df[profile_data[j].name], label=profile_data[j].fancy_name,color=plot_color)
		                    else:
		                        plt.plot(xval,profile_data[j].fac*df[profile_data[j].name], label=profile_data[j].fancy_name,color=plot_color)
		        title = Combined_plot_properties[sort-1].title+" after "+str(render_time[profile_sample[i]])+" "+t_unit[profile_sample[i]] 
		        plt.xlabel(xlab)
		        ylab = Combined_plot_properties[sort-1].ylab
		        plt.ylabel(ylab)
		        plt.suptitle(title)
		        ax.legend()
		        fn = Combined_plot_properties[sort-1].fname
		        filename = generic_title+fn+"_sample"+str(profile_sample[i])+".png"
		        fname = ''.join(filename)
		        print(fname)
		        #os.chdir(folder[f].fold)
		        plt.savefig(str(filename),dpi=500) #,transparent=True)
		        #os.chdir("..")

	#################### FLUX ##############################
	Flux_res = []
	flux_data_names = []

	for i in range(len(FluxVariables)):
		flux_data_names.append(FluxVariables[i].name)
	if(len(FluxVariables)>0):
		for i in range(len(folder)):
		    fold = folder[i].fold
		    Flux, fvarnames, nloc, fsTimes = extractFluxData(fold)
		    vec = np.fromstring(Flux[0],dtype=float, sep=" ")
		    numdata = np.zeros((len(Flux),len(vec)))
		    for i in range(len(Flux)):
		        vec = np.fromstring(Flux[i],dtype=float, sep=" ")
		        numdata[i][0] = 1+i % nloc
		        for j in range(1,len(vec)):
		            numdata[i][j] = vec[j]

		    df = pd.DataFrame(data=numdata,index=None,columns=fvarnames)
		    Flux_res.append(df)

		fplot_time = np.zeros(len(fsTimes))
		for t in range(len(fsTimes)):
		    fplot_time[t] = fsTimes[t]/3600/24/365


#	fnode_list = []
#	if(len(FluxVariables)>0):
		for i in range(len(fnode_list)):
		    for j in range(len(FluxVariables)):
		        fig, ax = plt.subplots()
		        for f in range(len(folder)):
		            dff = Flux_res[f]  
		            df=dff.loc[dff['Loc'] == fnode_list[i]]
		            if(len(folder)>1):
		                plot_color=folder[f].color
		            else:
		                plot_color=FluxVariables[j].color 
		            if(FluxVariables[j].acidity):
			            acidity = df["H[+]"]*0.001+df["HSO4[-]"]*0.097+2*df["H2SO4(aq)"]*0.098
			            plt.plot(fplot_time,acidity, label=str(folder[f].legend_label),color=plot_color)
		            else:
			            plt.plot(fplot_time,FluxVariables[j].fac*df[FluxVariables[j].name], label=str(folder[f].legend_label),color=plot_color)
		        
		        if(FluxVariables[j].settmax==True):
			        ax.set_xlim(right=fplot_time[FluxVariables[j].tmax])
		        if(FluxVariables[j].settmin==True):
			        ax.set_xlim(left=fplot_time[FluxVariables[j].tmin])
		        title = "Cumulative flux of "+FluxVariables[j].fancy_name+" at "+str(fnode_list[i])
		        xlab = "Time (years)"
		        plt.xlabel(xlab)
		        ylab = "Cumulative flux of "+FluxVariables[j].fancy_name+" ("+FluxVariables[j].unit+")" 
		        plt.ylabel(ylab)
		        plt.suptitle(title)
		        ax.legend()

		        filename = generic_title+"FluxCum_"+str(FluxVariables[j].fname)+".png"
		        fname = ''.join(filename)
		        #print(fname)
		        if(len(folder)==1):
		            #os.chdir(folder[f].fold)
		            plt.savefig(str(fname),dpi=500) #,transparent=True)
		            #plt.clf()
		            #os.chdir("..")
		        else:          
		            plt.savefig(str(fname),dpi=500) #,transparent=True)


class plot_var:

    logscale: bool
    fac: float
    acidity: bool
    setxmin: bool
    setymin: bool
    setxmax: bool
    setymax: bool
    settmin: bool
    settmax: bool
    tmin: float  
    tmax: float
    xmin: float
    ymin: float
    xmax: float
    ymax: float   # borner les axes xmin, xmax, ymin, ymax, tmin, tmax
    color: str
    name: str   # name: comme ça apparait dans hytec et les fichiers .res
    fancy_name: str  # pour renommer de manière "fancy"
    fname: str   # pour le nom du fichier qui va être généré
    transient: bool   # est-ce que tu veux que ce soit à un endroit au cour sdu temps
    profile: bool  # est-ce que veux profil spatial après un an
    to_plot:bool   # option inutile: 
    unit: str   # l'unité que tu peux spécifier de la grandeur (ex: mg/l)
    type_sort: int
        
    def __init__(self,name,unit="",fancy_name=" ",fname=" ",to_plot=True,logscale=False,fac=1,setxmin=False,acidity=False,setymin=False,setxmax=False, setymax=False,transient=True,type_sort=0,profile=True,xmin=-1,ymin=-1,xmax=-1, ymax=-1,settmin=False,settmax=False,tmin=-1,tmax=-1,color=""):
        self.logscale = logscale
        self.fac = fac
        self.setxmin = setxmin
        self.setymin = setymin
        self.setxmax = setxmax
        self.setymax = setymax
        self.settmin = settmin
        self.settmax = settmax
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax
        self.tmin = tmin
        self.tmax = tmax
        self.color = color
        self.type_sort = type_sort
        self.to_plot = to_plot
        if(self.xmin != -1):
            self.setxmin = True
        if(self.ymin != -1):
            self.setymin = True
        if(self.xmax != -1):
            self.setxmax = True
        if(self.ymax != -1):
            self.setymax = True
        if(self.tmax != -1):
            self.settmax = True
        if(self.tmin != -1):
            self.settmin = True
        self.acidity = acidity
        self.name = name
        self.unit = unit
        if(fancy_name == " "):
            self.fancy_name = self.name
        else:
            self.fancy_name = fancy_name
        if(fname == " "):
            self.fname = self.name
        else:
            self.fname = fname
        self.transient = transient
        self.profile = profile
    


class Simulation:
    fold: str
    legend_label: str
    plot_time: float
    color: str
    def __init__(self,fold,legend_label=" ", plot_time=[],color=""):
        self.fold = fold
        self.legend_label = legend_label
        self.plot_time = []
        self.color = color
        if(legend_label == " "):
            self.legend_label = self.fold
            
class Combined_plot:
    title: str
    logscale: bool
    ylab: str
    fname: str
    setxmin: bool
    setymin: bool
    setxmax: bool
    setymax: bool
    xmin: float
    ymin: float
    xmax: float
    ymax: float
    def __init__(self,title,fname,ylab,logscale=False,setxmin=False,setymin=False,setxmax=False, setymax=False,xmin=-1,ymin=-1,xmax=-1, ymax=-1):
        self.title = title
        self.fname = fname
        self.ylab = ylab
        self.logscale = logscale
        self.setxmin = setxmin
        self.setymin = setymin
        self.setxmax = setxmax
        self.setymax = setymax
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax
        if(self.xmin != -1):
            self.setxmin = True
        if(self.ymin != -1):
            self.setymin = True
        if(self.xmax != -1):
            self.setxmax = True
        if(self.ymax != -1):
            self.setymax = True

		

def extract_variable_name(line,doUnit):
    ind_begin = 0
    for j in range(len(line)):
        if(line[j] == "c" and line[j+1] == "o" and line[j+2] == "l" and line[j+3] == "u" and line[j+4] == "m" and line[j+5] == "n"):
            ind_begin = j+10
            if(line[j+10] == ":"):
                ind_begin = j+11
    if(ind_begin == 0):
        varstr = "0"
        unit = "0"
        return varstr, unit
    
    for j in range(ind_begin+1,len(line)):
        if(line[j].isspace()):
            ind_end = j
            break
    if(doUnit):
        for j in range(ind_end+1,len(line)):
            if(line[j].isspace() == False):
                unit_begin = j
                break
            j+=1
        unit_n = line[unit_begin:len(line)-1]
        unit = ''.join(unit_n)
    else:
        unit = "0"
    var = line[ind_begin:ind_end]
    
    if(var[0].isspace()):
        varr = var[1:len(var)]
    else:
        varr = var
    varstr = ''.join(varr)
    return varstr, unit
    


def extract_sample_time(line):
    ind_begin = 0
    for j in range(len(line)):
        if(line[j] == "t" and j<len(line)-6):
            if(line[j] == "t" and line[j+1] == "i" and line[j+2] == "m" and line[j+3] == "e"):
                ind_begin = j+7
    if(ind_begin == 0 and line[0].isspace() == True):
        time = -1
        return time
    elif(line[2] == "s" and line[3] == "e" and line[4] == "g"):
        time = -1
        return time
    elif(ind_begin == 0):
        time = -2
        return time
    for j in range(ind_begin+1,len(line)):
        if(line[j].isspace() or line[j] == "\t"):
            ind_end = j
            break
    var = line[ind_begin:ind_end]
    time = float(var)
    return time

def extend_vector(a,n,l,m):
	vec = np.zeros((n,l))
	for i in range(l):
		for j in range(m):
			vec[i][j] = a[i][j]
	return vec

def extractData(fold):
	hytec_file = os.path.join(fold, "HYTEC.res")
	hytec_data = pd.read_csv(hytec_file)
	i = 0
	inPreamble = True
	varNames = []
	unitNames = []
	nvar = 0
	while(inPreamble):
		line = hytec_data.iloc[[i]]
		str_line = line.to_string(header=False) 
		newvar, unit = extract_variable_name(str_line,True)
		if(newvar == "0" and nvar > 0):
			inPreamble = False
		elif(newvar != "0"):
			varNames.append(newvar)
			unitNames.append(unit)
			nvar+=1        
		i=i+1 
	#varNames.append("time_s")
	sampleTimes = []
	nsample = 0
	nosample = False
	first = False
	while(i<hytec_data.shape[0]):
		line = hytec_data.iloc[[i]]
		str_line = line.to_string(header=False) 
		newvar = extract_sample_time(str_line)
		if(newvar != -1 and nosample == False):
			if(first==False):
				#first = True
				first_data_line = i+1
			i_begin = i+1
			nosample = True
			sampleTimes.append(newvar)
			nsample+=1
		elif((newvar >= 0 and nosample == True) or i==hytec_data.shape[0]-1):
			if(first==False):
				first = True
				last_data_line = i-2
			i_end = i-2
			if(i==hytec_data.shape[0]-1):
				i_end = i
			else:
				i = i - 1
			nosample = False        
		i+=1
	Nnodes = last_data_line-first_data_line+1
	#os.chdir("..")
	return varNames, sampleTimes, Nnodes, first_data_line, unitNames


def get_profile_data(datalist,samples,fold,firstLine,Nn,varNames):
	hytec_file = os.path.join(fold, "HYTEC.res")
	hytec_data = pd.read_csv(hytec_file)
	nvar = len(varNames)
	T = samples
	begin = firstLine+T*(Nn+2)
	end = begin+Nn
	data = np.zeros((end-begin,nvar))  
	nl = 0
	for j in range(begin,end):
		numline = hytec_data.iloc[[j]].to_string(header=False)
		vec = np.fromstring(numline,dtype=float, sep=" ")
		data[nl][0:len(vec)-1] = vec[1:len(vec)]
		nl += 1
	newdf = pd.DataFrame(data=data,index = None, columns = varNames)
	shortdf = newdf[datalist]

	return shortdf


def get_transient_data(datalist,fold,firstLine,node,Nn,varNames,times):
	Nsamples = len(times)
	hytec_file = os.path.join(fold, "HYTEC.res")
	hytec_data = pd.read_csv(hytec_file)
	nvar = len(varNames)
	#datalist.append("time")
	#varNames.append("time")
	data = np.zeros((Nsamples,nvar))
	for i in range(Nsamples):
		numline = hytec_data.iloc[[firstLine+node-1+i*(Nn+2)]].to_string(header=False)
		vec = np.fromstring(numline,dtype=float, sep=" ")
		data[i][0:len(vec)-1] = vec[1:len(vec)]
		#data[i][len(vec)-1]=times[i]
		i+=1        
	newdf = pd.DataFrame(data=data,index = None, columns = varNames)
	shortdf = newdf[datalist]

	return shortdf
                






def extractFluxData(fold):
	hytec_file = os.path.join(fold, "HYTEC.res")
	hytec_data = pd.read_csv(hytec_file)
	i = 0
	inPreamble = True
	varNames = ["Loc"]
	nvar = 0
	while(inPreamble):
		line = hytec_data.iloc[[i]]
		str_line = line.to_string(header=False) 
		newvar, useless_units = extract_variable_name(str_line,False)
		if(newvar == "0" and nvar > 0):
			inPreamble = False
		elif(newvar != "0"):
			varNames.append(newvar)
			nvar+=1        
		i=i+1 
	#varNames.append("time_s")
	sampleTimes = []
	nsample = 0
	nosample = False
	first = False
	SampleList = []
	bigDf = []
	newdf = pd.DataFrame(index = None, columns = varNames)
	while(i<hytec_data.shape[0]):
		line = hytec_data.iloc[[i]]
		str_line = line.to_string(header=False) 
		newvar = extract_sample_time(str_line)
		if(newvar >= 0 and nosample == False):            
			i_begin = i+1
			nosample = True
			sampleTimes.append(newvar)
			nsample+=1
		elif((newvar >= 0 and nosample == True) or i==hytec_data.shape[0]-1):
			i_end = i-2
			if(i==hytec_data.shape[0]-1):
				i_end = i
				nloc = i_end-i_begin+1
			else:
				i = i - 1
			nosample = False
			#data = np.zeros((i_end-i_begin+1,nvar))
			#nl = 0
			for j in range(i_begin,i_end+1):
				numline = hytec_data.iloc[[j]].to_string(header=False)
				l_numline = list(numline)
				stringline = ""
				for k in range(len(l_numline)): 
					if(l_numline[k] != "\t" and l_numline[k] != "t" and l_numline[k] != " "):
						stringline+=l_numline[k]
				laststring = ""
				for k in range(len(stringline)):
					if(stringline[k].isalnum()==False and stringline[k]!="." and stringline[k]!="-"):                        
						laststring += ' '
					else:
						laststring += stringline[k]
				numline = str(laststring)
				vec = np.fromstring(numline,dtype=float, sep=" ")
				bigDf.append(numline)
				
		i+=1
	#newdf.head()
	#newdf = pd.DataFrame(data=bigDf,index = None, columns=varNames, sep = " ")
	return bigDf, varNames, nloc, sampleTimes

if __name__== '__main__':
	main()






