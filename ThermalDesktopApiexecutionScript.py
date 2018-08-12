#Thermal Desktop API execution Script


import sys 
import sys.path

import clr
TAV = clr.Addreference('TdApiV1')#This is the imported TdApiV1 module

#Create Instance of Thermal Desktop
td = TdApiV1.ThermalDesktop() #instantiates an instance of thermal desktop object
td.Connect() # this turns on autodesk and thermal desktop
td.Print("Hello World") #This prints Hello World to the Autodesk Command Interface

# Add Node Properties

submod = td.CreateSubmodel('HelloWorld')#creates a submodel named HelloWorld
node1.Submodel = submod.get_Name()
node2.Submodel = submod.get_Name()

#Create a Connection Object For Both objects that will be connected
con1 = TdApiV1.Connection()
con1.set_Handle(node1.Handle)
con2 = TdApiV1.Connection()
con2.set_Handle(node2.Handle)


cond1 = td.CreateConductor(con1.Handle, con2.Handle)#cond numbering is arbitrary
cond1.Submodel = submod.get_Name()#Set submodel name
cond1.ValueExp.set_Value("100.0 + 100.0")


Units = TdApiV1.Units
cond1.ValueExp.set_units(TdApiV1.UnitsData(Units.SI))


td.SetRcConductor(con1)


td.ZoomExtents()



