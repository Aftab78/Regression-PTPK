'''
Created on November 18, 2022

@author: AFTAB HASSAN
'''

import csv,os
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from copy import deepcopy
class PtpkAnalysis:
    def __init__(self):        
        
        self.plants={} # {('KRI','Rail):'Krishnapatnam,('KRI','Road):'Krishnapatnam}
        self.allData=[]
        self.rsquaredResult=[]
        self.output=[]   
        self.Path="OutputImages"
        self.slabType=["Slab-1","Slab-2","Slab-3","Slab-4","Slab-5"] # Slab-1 : 0-250 ,Slab-2 : 250-500 ,Slab-3 : 500-750 ,Slab-4 : 750- up to, else Slab-5
        self.slabRange=[0,800,1500]
        self.intervalValue=100
        self.SlabPercentage=[50,30,20]
        self.slabToggle=0  # 0 -> OFF & 1 -> ON
        self.rangeToggle=0 # 0 -> slabRange according to Distance & 1 -> SlabPercentage according to frequency percentage
        
        try:
            os.makedirs(self.Path, exist_ok = True)
            print("Directory '%s' created successfully" % self.Path)                      
        except OSError as error:
            print("Directory '%s' can not be created" % self.Path) 
            print("error" % error)
         
        self.loadInputData()  
         
    def objective(self,x, a, b):
            return a * x ** b
        
    def intervalFrequency(self,distanceList,intervalValue):
        frequency = {}
        totdata=len(distanceList)
        fromValue=0
        toValue=intervalValue
        numberOfSlab=int(totdata)
          
        loopCount=0
        while(len(distanceList)!=0):
            temp=[]
            
            if(loopCount<(numberOfSlab-1)):
                #print("In if")
                for d in distanceList:
                    if (d<=toValue):
                        temp.append(d)
                        
                for d in temp:
                    distanceList.pop(distanceList.index(d))
                            
                frequency[(fromValue,toValue)]=[len(temp),round(len(temp)/totdata*100,2)]
                fromValue=toValue
                toValue+=intervalValue
                loopCount+=1
                
            else: 
                #print("In else")        
                for d in distanceList:        
                    temp.append(d)            
                    
                frequency[(fromValue,max(temp))]=[len(temp),round(len(temp)/totdata*100,2)]
                fromValue=toValue
                toValue=max(temp)  
                loopCount+=1   
                for d in temp:
                    distanceList.pop(distanceList.index(d))            
              
        frequencyTemp=deepcopy(frequency)
        slabRanges=[]
        minInterval,maxInterval=list(frequencyTemp.keys())[0]
        
        for slab in range(len(self.SlabPercentage)):
            percent=0
            temp=[]
            for interval,countPercentage in frequencyTemp.items():
                #print(interval[0],"\t",interval[1],"\t",countPercentage[0],"\t",countPercentage[1])
                percent+=countPercentage[1]
                temp.append(interval)
                maxInterval=interval[1]
                if(percent>self.SlabPercentage[slab]):
                    #print("Percentage : ",percent,minInterval,maxInterval)
                    slabRanges.append(minInterval)
                    minInterval=maxInterval
                    break
            if(slab==len(self.SlabPercentage)-1):
                #print("Percentage : ",percent,minInterval,maxInterval)
                slabRanges.append(minInterval)
                minInterval=maxInterval
                    
            for key in temp:
                frequencyTemp.pop(key)
        
            
        return slabRanges
        
    def slabBuilder(self,fileName):
        #print("In Slab range builder.")  
        file = open(fileName)
        csvdata = csv.reader(file)
        next(csvdata) # for skip Header
        distanceList=[]
        for row in csvdata:                
            distanceList.append(float(row[8]))
            
        #print(distanceList)      
        if(self.rangeToggle==1):
            self.slabRange=self.intervalFrequency(distanceList,self.intervalValue)        
        
        
        
    def loadInputData(self):
        #print("In read_csv_data ")
        fileName='data.csv'
        file = open(fileName)
        csvreader = csv.reader(file)
        header = next(csvreader)
        print(header)
        
        if(self.slabToggle==1):
            self.slabBuilder(fileName)            
        
        for row in csvreader:
            if(self.slabToggle==1):
                
                distance=float(row[8])
                #print("distance = ",distance)
                # *** Slab Type , for Plant according to Distance
                for indx in range(len(self.slabRange)):
                    if(indx<(len(self.slabRange)-1)):
                        if(distance>self.slabRange[indx] and distance<=self.slabRange[indx+1]):
                            slab=self.slabType[indx]
                    else:
                        # Last Slab
                        if(distance>self.slabRange[indx]): 
                            slab=self.slabType[indx]                
                    
                    
                #print(slab)   
                row.append(slab)
                self.plants[(row[0],row[4],row[5],row[6],row[9])]=row[1]
                self.allData.append(row)
                
            else:
                self.plants[(row[0],row[4],row[5],row[6])]=row[1]
                self.allData.append(row)
          
    def plotGraph(self,a,b,eq,plotX,plotY,R_square,plantName,mode,slab,maxy_limit):
        # Create the vectors X and Y
        x = np.array(range(1,3000))
        y = self.objective(x, a, b)
        #plt.ylim([0,50])
        plt.ylim([0,maxy_limit])
        
        # Create the plot
        if(len(plotX)>=4):
            plt.plot(x,y,label='y='+eq+'\n R^2 = '+str(R_square))
        else:
            plt.plot(x,y,label='y='+'In sufficient data'+'\n R^2 = '+str(R_square))
        
        # plotting points as a scatter plot
        plt.scatter(plotX, plotY, color= "green",marker= ".", s=100)
        
        # Add a title
        plt.title('PTPK Graph - '+plantName+' - '+mode+' - '+slab)
        
        # Add X and y Label
        plt.xlabel('Distance')
        plt.ylabel('PTPK ( pmt / distance )')
        
        # Add a grid
        plt.grid(alpha=.4,linestyle='--')        
        # Add a Legend
        plt.legend()
        # Show the plot
        plt.savefig(self.Path+'\\PTPK_Graph_'+plantName+'-'+mode+'-'+slab+'.jpg',facecolor = '#1589FF', orientation='landscape',dpi=600)
        plt.cla()
        #plt.show()
        

    def output_csv(self):        
        
        fp=open("RsquareResult"+".csv","w")
        
        header='Plant Name'+','+'Mode'+','+'Slab Type'+','+'Equation'+','+'RSquare'                
        fp.writelines(header)
        fp.write('\n') 
        
        for x in self.rsquaredResult:  
            #print(x)         
            temp=""         
            for i in range(len(x)-1):
                temp=temp+str(x[i])+','
                
            temp=temp+str(x[len(x)-1])        
            fp.writelines(temp)
            fp.write('\n')

        fp.close()
        
        fp1=open("FreightResult"+".csv","w")
        
        header='Plant Name'+','+'Mode'+','+'Destination'+','+'PMT'+','+'Distance'+','+'PTPK'+','+'Derived'+','+'Deviation'+','+'Freight'+','+'Slab Type'       
        fp1.writelines(header)
        fp1.write('\n')        
       
        for x in self.output:           
            temp=""         
            for i in range(len(x)-1):
                temp=temp+str(x[i])+','
                
            temp=temp+str(x[len(x)-1])        
            fp1.writelines(temp)
            fp1.write('\n')

        fp1.close()

        print ('File generated successfully.')
    
    def RegressionAnalaysisLogic(self):
                    
        #print("Plant and Mode : ",self.plants)
        for plantModeCode,plantName in self.plants.items():
            
            PlantCode=plantModeCode[0]           
            mode=plantModeCode[1]
            LinkType=plantModeCode[2]
            RouteType=plantModeCode[3]
            slabType="Slab-1"
            #print(PlantCode,mode,plantName)
            plantData=[]  
            plotX=[]
            plotY=[]
            actual=[]
            predict=[]
            coeff_a=1   ## a in equation
            coeff_b=1   ## b in equation            
            
            for row in self.allData:
                if(self.slabToggle==1):
                    slabType=plantModeCode[4]
                    if(row[0]==PlantCode and row[4]==mode and row[5]==LinkType and row[6]==RouteType and row[9]==slabType):
                        plantData.append([row[2],float(row[7]),float(row[8])])
                        plotX.append(float(row[8])) # Distance
                        plotY.append(float(row[7])/float(row[8])) # PTPK = pmt/distance
                else:
                    if(row[0]==PlantCode and row[4]==mode and row[5]==LinkType and row[6]==RouteType):
                        plantData.append([row[2],float(row[7]),float(row[8])])
                        plotX.append(float(row[8])) # Distance
                        plotY.append(float(row[7])/float(row[8])) # PTPK = pmt/distance
                        
            #print(plantData)
            if(len(plotX)>=4):
                generatedEquation,coeff=self.generate_equation(plotX,plotY) 
                print("Generated Equation : y=",generatedEquation)
                  
                coeff_a=float(coeff[0])
                coeff_b=float(coeff[1]) 
            else:
                generatedEquation="In sufficient data & data count should be >=4"
                print("Generated Equation : ",generatedEquation)             
            
            maxy_limit=-99999            
            for data in plantData:                
                destination=data[0]
                pmt=data[1]
                distance=float(data[2])
                ptpk=pmt/distance
                if(len(plotX)>=4):
                    derivedPtpk=self.y_equation(coeff_a,distance,coeff_b)
                else:
                    derivedPtpk=pmt/distance
                deviation=abs(ptpk-derivedPtpk)
                freightPmt=derivedPtpk*distance
                if(maxy_limit<ptpk):
                    maxy_limit=ptpk                
                    
                #print(plantName,mode,destination,pmt,distance,ptpk,derivedPtpk,deviation,freightPmt)
                self.output.append([plantName,mode,destination,pmt,distance,ptpk,derivedPtpk,deviation,freightPmt,slabType])
        
                actual.append(ptpk)
                predict.append(derivedPtpk)
            
            R_square=round(self.RSquared(actual,predict),5)
            print("R^2 = ",R_square)     
            self.rsquaredResult.append([plantName,mode,slabType,'y='+generatedEquation,R_square])
            self.plotGraph(coeff_a,coeff_b,generatedEquation,plotX,plotY,R_square,plantName,mode,slabType,maxy_limit+2)
            
        self.output_csv()
        
    def y_equation(self,a,x,b):              
        return self.objective(x,a,b)
    
    def RSquared(self,actual,predict):
        corr_matrix = np.corrcoef(actual, predict)
        corr = corr_matrix[0,1]
        R_sq = corr**2      
        return R_sq
    
    def generate_equation(self,plotx,ploty):
               
        # load the dataset        
        x, y = plotx,ploty
        # curve fit
        popt, _ = curve_fit(self.objective, x, y)     
        a, b = popt
        equ='%.3fx^%.3f' % (a, b)
        #print('y = %.3fx^%.3f' % (a, b))
        
        return equ,(a,b)     

def main():    
    Ptpk = PtpkAnalysis()
    Ptpk.RegressionAnalaysisLogic()
    
if __name__ == '__main__':
    main()