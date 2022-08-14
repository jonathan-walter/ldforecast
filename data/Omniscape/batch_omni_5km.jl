#Batch-run omniscape using multiple .ini files

#Defaults to using 1 thread. For parallel processing, see
#Extensions -> Julia -> settings (gear icon) and scroll down to Num threads

using Omniscape #launches omniscape

#sets working directory -- REMEMBER TO CHANGE IF NEEDED
cd("/Users/kathryngrage/Box Sync/Gypsy Moth/Jeffress Spread Forecasting/JuliaScripts")

#run_omniscape("./testing/test.ini") #runs test 
#run_omniscape("./testing/test2.ini") #runs test2

run_omniscape("./Full Extent/2017_5km.ini") 
run_omniscape("./Full Extent/2018_5km.ini") 
run_omniscape("./Full Extent/2019_5km.ini") 
#run_omniscape("./Full Extent/2016_5km.ini") 
#run_omniscape("./Full Extent/2015_5km.ini") 
#run_omniscape("./Full Extent/2014_5km.ini") 
#run_omniscape("./Full Extent/2013_5km.ini") 
#run_omniscape("./Full Extent/2012_5km.ini") 
#run_omniscape("./Full Extent/2011_5km.ini") 
#run_omniscape("./Full Extent/2010_5km.ini") 
#run_omniscape("./Full Extent/2009_5km.ini") 
#run_omniscape("./Full Extent/2008_5km.ini") 
#run_omniscape("./Full Extent/2007_5km.ini") 
#run_omniscape("./Full Extent/2006_5km.ini") 
#run_omniscape("./Full Extent/2005_5km.ini") 
#run_omniscape("./Full Extent/2004_5km.ini")
#run_omniscape("./Full Extent/2003_5km.ini") 
#run_omniscape("./Full Extent/2002_5km.ini") 
#run_omniscape("./Full Extent/2001_5km.ini") 