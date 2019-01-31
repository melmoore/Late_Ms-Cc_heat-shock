#Ms+Cc Late larval heat shock--PRELIMINARY PLOTS AND ANALYSIS

#----------------

#load libraries

library(readr)
library(ggplot2)
library(tidyr)
library(Rmisc)
library(dplyr)
library(viridis)
library(cowplot)
library(extrafont)


#---------

#Load data


lhs <- read_csv("Ms+Cc_late_hs_incomp_8-10-18.csv", 
                col_types = cols(hs.treat = col_factor(levels = c("0", "hs.4", "hs.5")), 
                                 temp.var = col_factor(levels = c("0", "10"))))
View(lhs)


#-------------

#remove empty rows and columns

lhs$id[is.na(lhs$id)]<-0

lhs<-subset(lhs,id!=0)

lhs<-lhs[,-(100:112)]


#-----------------


#remove dead caterpillars

lhs$date.died[is.na(lhs$date.died)]<-0

lhs.cl<-subset(lhs, date.died=="0")


#--------------------

#Making a long data set of mass48em and mass at end of feeding trial

feed.lng<-gather(lhs.cl, time, mass, mass.48em, mass.postem.end)
View(feed.lng)

#subset out cats with very large masses (>10,000mg), probably errors. fix in cleaning

feed.lng<-subset(feed.lng,mass<10000)


#plotting mass at em and mass at end of feeding trial

postem.mass.plot<-ggplot(feed.lng, aes(x=time, y=mass, group=interaction(hs.treat,id), color=hs.treat))
postem.mass.plot+geom_point(
               )+geom_line(
               )+facet_wrap(~temp.var)



#Plotting mean mass for feeding trial

feed.sum<-summarySE(feed.lng, measurevar = "mass",
                    groupvars = c("temp.var", "hs.treat","time"),
                    na.rm = TRUE)
feed.sum


mnfeed.plot<-ggplot(feed.sum, aes(x=time, y=mass, group=hs.treat, color=hs.treat))
mnfeed.plot+geom_point(
)+geom_line(
)+facet_wrap(~temp.var)





#-------------

#plotting female sex ratio (just of ecl adults, and not including late ecl)

lhs.cl$fem.rat<-lhs.cl$live.fem/lhs.cl$num.ecl

fem.sum<-summarySE(lhs.cl, measurevar = "fem.rat",
                   groupvars = c("temp.var", "hs.treat"),
                   na.rm = TRUE)
fem.sum


fem.plot<-ggplot(lhs.cl, aes(x=hs.treat, y=fem.rat, group=temp.var, color=temp.var))
fem.plot+geom_jitter()

mnfem.plot<-ggplot(fem.sum, aes(x=hs.treat, y=fem.rat, group=temp.var, color=temp.var))
mnfem.plot+geom_point(
)+geom_errorbar(aes(ymin=fem.rat-se, ymax=fem.rat+se))



#------------------------------

#Calculating development time

#converting dates to julian dates

##Converts x into julian date
j.date<-function(x){
  strptime(x, "%m/%d")$yday+1
}

#Takes all columns that have "date." in the name, and converts contents to Julian day using j.date function. Renames columns (adds a
##j to end of column name), and binds the out put julian day columns to the original data set

lapj.date<-function(df){
  date.j<-lapply(df[,grep("date.",colnames(df))],j.date)
  date.j<-as.data.frame(date.j)
  colnames(date.j)<-paste(colnames(date.j), "j", sep = ".")
  output.df<-cbind(df,date.j)
  output.df
}


lhs.cl<-lapj.date(lhs.cl)



#Calculating ages

lhs.cl$age.3<-lhs.cl$date.3.j-lhs.cl$date.hatch.j
lhs.cl$age.4<-lhs.cl$date.4.j-lhs.cl$date.hatch.j
lhs.cl$age.5<-lhs.cl$date.5.j-lhs.cl$date.hatch.j
lhs.cl$age.48em<-lhs.cl$date.em.j-lhs.cl$date.hatch.j
lhs.cl$age.postem.end<-lhs.cl$date.postem.7.j-lhs.cl$date.hatch.j

#-------------------------

#Plotting mass X age for hosts

#Making long dataset of mass

lhs.lng<-gather(lhs.cl, instar, mass, mass.3, mass.4, mass.5, mass.48em, mass.postem.end)
lhs.lng$instar<-gsub("mass.", "", lhs.lng$instar)

#Making long dataset of age

lhs.age<-gather(lhs.cl, instar, age, age.3, age.4, age.5, age.48em, age.postem.end)
lhs.age$instar<-gsub("age.", "", lhs.age$instar)


#Remove extra columns from lhs.age

lhs.age<-select(lhs.age, id, hs.treat, instar, age)


#merge age and mass data frames

lhs.lng<-merge(lhs.lng, lhs.age, by=c("id", "hs.treat", "instar"))


#Toss individuals with mass over 10,000--probably typos

lhs.lng<-subset(lhs.lng, mass<10000)

#Finding log mass

lhs.lng$mass<-as.numeric(lhs.lng$mass)

lhs.lng$log.mass<-log(lhs.lng$mass)



#----------------------

#plotting mass X age

mage.plot<-ggplot(lhs.lng, aes(x=age, y=log.mass, group=interaction(id, hs.treat), color=hs.treat))
mage.plot+geom_point(
)+geom_line(
)+facet_wrap(~temp.var)



#Finding summarySE of age and mass

lm.sum<-summarySE(lhs.lng, measurevar = "log.mass",
                  groupvars = c("temp.var", "hs.treat", "instar"),
                  na.rm=TRUE)
lm.sum


age.sum<-summarySE(lhs.lng, measurevar = "age",
                   groupvars = c("temp.var", "hs.treat", "instar"),
                   na.rm=TRUE)
age.sum


lm.sum$age<-age.sum[,5]
lm.sum$age.se<-age.sum[,7]


mage.sum.plot<-ggplot(lm.sum, aes(x=age, y=log.mass, group=hs.treat, color=hs.treat))
mage.sum.plot+geom_point(size=3
)+geom_line(size=1.2
)+geom_errorbar(aes(ymin=log.mass-se, ymax=log.mass+se)
)+facet_wrap(~temp.var)


mage.sum.plot2<-ggplot(lm.sum, aes(x=age, y=log.mass, group=temp.var, color=temp.var))
mage.sum.plot2+geom_point(size=3
)+geom_line(size=1.2
)+facet_wrap(~hs.treat)


#---------------------------------

#plotting wasp development time

lhs.cl$ttem<-lhs.cl$date.em.j-lhs.cl$date.ovp.j
lhs.cl$ttecl<-lhs.cl$date.ecl.j-lhs.cl$date.ovp.j

devem.sum<-summarySE(lhs.cl, measurevar = "ttem",
                     groupvars = c("temp.var", "hs.treat"),
                     na.rm = TRUE)
devem.sum


devecl.sum<-summarySE(lhs.cl, measurevar = "ttecl",
                      groupvars = c("temp.var", "hs.treat"),
                      na.rm = TRUE)
devecl.sum


devem.sum$ttecl<-devecl.sum[,4]
devem.sum$ttecl.se<-devecl.sum[,6]


ccdevem.plot<-ggplot(devem.sum, aes(x=temp.var, y=ttem, group=hs.treat, color=hs.treat))
ccdevem.plot+geom_point(size=3
)+geom_line(size=1.2)


ccdevecl.plot<-ggplot(devem.sum, aes(x=temp.var, y=ttecl, group=hs.treat, color=hs.treat))
ccdevecl.plot+geom_point(size=3
)+geom_line(size=1.2
)+geom_errorbar(aes(ymin=ttecl-se, ymax=ttecl+se))



#------------------------

#Plotting the number of eclosed wasps

numecl.plot<-ggplot(lhs.cl, aes(x=hs.treat, y=tot.ecl, group=temp.var, color=temp.var))
numecl.plot+geom_jitter(size=3)


numecl.sum<-summarySE(lhs.cl, measurevar = "tot.ecl",
                      groupvars = c("temp.var", "hs.treat"),
                      na.rm = TRUE)
numecl.sum

numeclsum.plot<-ggplot(numecl.sum, aes(x=temp.var, y=tot.ecl, group=hs.treat, color=hs.treat))
numeclsum.plot+geom_point(size=3
)+geom_line(size=1.2
)+geom_errorbar(aes(ymin=tot.ecl-se, ymax=tot.ecl+se), width=.5)



#----------------------------

#looking at number em and mass at em

#exclude super large mass48em number--probs typo

lhs.cl<-subset(lhs.cl, mass.48em<10000)

massem.plot<-ggplot(lhs.cl, aes(x=mass.48em, y=num.em, group=hs.treat, color=hs.treat))
massem.plot+geom_point(
)+geom_smooth(method=lm, se=TRUE
)+facet_wrap(~temp.var)


#------------------------------

#Making long data frame of wasp pupal stages

stage.lng<-gather(lhs.cl, stage, num.stage, live.fem, live.male, dead.fem, dead.male, dead.larv, dead.epup, dead.mpup.tot, dead.lpup.tot)


#Plotting the number in each stage for hs.treat and temp.var

stage.plot<-ggplot(stage.lng, aes(x=tot.ecl, y=num.stage, group=hs.treat))
stage.plot+geom_point(
)+geom_smooth(method=lm
)+facet_wrap(~stage*temp.var)


#trying boxplots of stages

stbx.plot<-ggplot(stage.lng, aes(x=stage, y=num.stage))
stbx.plot+geom_boxplot(aes(color=hs.treat)
)+facet_wrap(~temp.var)


#Removing live individuals and outlier with very high num.stage

st.dead.lng<-subset(stage.lng, stage!="live.fem")
st.dead.lng<-subset(st.dead.lng, stage!="live.male")
st.dead.lng<-subset(st.dead.lng, num.stage<90)

#reordering factor levels
x = factor(x,levels(x)[c(4,5,1:3)])

st.dead.lng$stage<-factor(st.dead.lng$stage, levels=c("dead.larv", "dead.epup", "dead.mpup.tot", "dead.lpup.tot", "dead.fem", "dead.male"))


#trying boxplots of stages with only dead wasps

stbx.plot2<-ggplot(st.dead.lng, aes(x=stage, y=num.stage))
stbx.plot2+geom_boxplot(aes(fill=temp.var)
)+facet_wrap(~hs.treat)


#Trying as percentage of tot.dead

st.dead.lng$perc.stage<-st.dead.lng$num.stage/st.dead.lng$tot.dead

stbx.plot2<-ggplot(st.dead.lng, aes(x=stage, y=perc.stage))
stbx.plot2+geom_boxplot(aes(fill=temp.var)
)+facet_wrap(~hs.treat)



stbx.plot3<-ggplot(st.dead.lng, aes(x=stage, y=perc.stage))
stbx.plot3+geom_boxplot(aes(fill=hs.treat)
)+facet_wrap(~temp.var)





#-------------
#plotting perc.surv for those that have been dissected

lhs.dis<-subset(lhs.cl, num.unem!="0")

lhs.dis$perc.surv<-lhs.dis$tot.ecl/lhs.dis$load


percsurv.plot<-ggplot(lhs.dis, aes(x=load, y=perc.surv, group=hs.treat, color=hs.treat))
percsurv.plot+geom_point(
)+geom_smooth(method=lm
)+facet_wrap(~temp.var)

