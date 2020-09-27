```R
# Fish.automatic
Automatic generation of fishplots from an event matrix.
#######################################
#Automatic fishplot generation from EM#
#######################################
# install.packages("devtools")
# install_github("chrisamiller/fishplot")

library(devtools)
library(fishplot)

#Load the data set.
load_matrix <- function(filename, sheetname) {
  data <- as.data.frame(read_xlsx(filename, sheetname))
  subdata <- data[ c(1:nrow(data)), c(1:ncol(data)) ]
  subdata <- as.matrix(subdata[complete.cases(subdata),])
  return(subdata)
}

EM <- load_matrix(filename="EM.xlsx",sheetname ="Sheet3") #Loading the EM.

#Assign which cells belong to the same time point.
samples <- list(2:13, #O
                49:59,49:59, #G1
                60:63,60:63,60:63, #G7
                28:48,28:48, #lp
                14:27)#hp

cells_tot <- list() #Making a list for the total cell number.

i <- 1
for(i in 1:length(samples)){

  columns <-unlist(samples[i])
  cells_tot[[i]] <- sum(as.numeric(EM[2,columns]))

  i <- i+1
}

#Provide a list of timepoints to plot.
#You may need to add fabricated points to end up with the desired visualization.
timepoints <- c(0,      #Start
                100,150, #G1
                200,    #Bottleneck after G1.
                250,300, #G7
                350,    #Bottleneck after G7.
                400,600) #lp and hp.

#Calculate each event's prevalence at each time point.
i <- 3
sum_row <- matrix(0,nrow(EM),length(samples)*2+1)
for(i in 3:nrow(EM)){

  j <- 1
  for(j in 1:length(samples)){
    print(j)
    samples_unl <- unlist(samples[j])
    themax <- max(samples_unl)
    themin <- min(samples_unl)

    sum_row[i,1] <- EM[i,1]
    sum_row[i,j+j] <- sum(as.numeric(EM[i,themin:themax])*as.numeric(EM[2,themin:themax]))
    sum_row[i,j+j+1] <- as.numeric(sum_row[i,j+j])/as.numeric(unlist(cells_tot[j]))*100

  }

  i <- i+1
}

#Making the final fract.table by extracting the fractions from sum_row.
#The odd columns contain the fraction of cells harboring each alteration.
#The even columns contain the number of cells harboring each alteration.
i <- 1
for(i in 1:(length(timepoints)+1)){
  if(i == 1){
    frac.table_original <- sum_row[3:nrow(sum_row),i+i-1]
  }else{
    frac.table_original <- cbind(frac.table_original,sum_row[3:nrow(sum_row),i+i-1])
  }
}

#Ordering the matrix a bit.
ordering <- frac.table_original[,2:ncol(frac.table_original)]
class(ordering) <- "numeric"
frac.table_original_ordering <- cbind(rowSums(ordering),frac.table_original)
frac.table_original <- frac.table_original_ordering[order(-as.numeric(frac.table_original_ordering[,1])),]
frac.table_original <- frac.table_original[frac.table_original[,1]>3,]
frac.table_original <- frac.table_original[,2:ncol(frac.table_original)]

parents <- matrix(0,ncol(frac.table_original),nrow(frac.table_original))
parents[1,] <- frac.table_original[,1]

#Determining the parental order of events.
t <- 2
for(t in 2:(length(timepoints)+1)){
  print("Här är t")
  print(t)
  frac.table <- frac.table_original[order(-as.numeric(frac.table_original[,t])),]

  i <- 1
  for(i in 1:nrow(frac.table)){

    if(i == 1){
      parent_pos <- match(frac.table[1,1],parents[1,])
      parents[t,parent_pos] <- 0
    }else{
      if(frac.table[i,t]!="0"){

        name_d <- frac.table[i,1] #Name of the daughter.
        row_d <- match(name_d,EM[,1]) #The row of the daughter in the EM.
        cells_d <- as.matrix(EM[,EM[row_d,]=="1"]) #Finding all columns having this daughter event.
        #cells_d_tot <- sum(as.numeric(cells_d[row_d,1:ncol(cells_d)]))

        j <- 1
        mother_found <- "0"
        while(mother_found == "0" && j < i){

          name_m <- frac.table[i-j,1] #Name of the potential mother.
          row_m <- match(name_m,EM[,1]) #The row of the mother in the EM.
          #cells_m_tot <- sum(as.numeric(cells_d[i-j,1:ncol(cells_d)]))

          difference <- as.numeric(frac.table[i-j,2:ncol(frac.table)])-as.numeric(frac.table[i,2:ncol(frac.table)])
          difference <- (difference >= 0)
          print(difference)
          if(all(difference,TRUE) == TRUE){
            #It is larger than the daughter in all samples.

            prevalent <- as.numeric(EM[row_m,2:ncol(EM)])-as.numeric(EM[row_d,2:ncol(EM)])

            prevalent <- (prevalent >= 0)
            print(prevalent)

            if(all(prevalent,TRUE)==TRUE){
              #The mother exist in all samples that the daughter does.
              daughter_col <- match(name_d,parents[1,])
              mother_col <- match(name_m,parents[1,])
              parents[t,daughter_col] <- mother_col
              mother_found <- "1"

            }

          }
          j <- j+1
        }

      }else{parents[t,match(frac.table[i,1],parents[1,])] <- "NP"}
    }
    i <- i+1
  }
  t <- t+1
}

frac.table_new <- frac.table_original[,2:ncol(frac.table_original)]
class(frac.table_new) <- "numeric"
frac.table_new <- round(frac.table_new,0)

parents_new <- matrix(0,1,ncol(parents))

#Determining the common parent between all timepoints.
i <- 1
for(i in 1:ncol(parents)){
  #mothers <- as.vector(as.numeric(parents[2:nrow(parents),i]))
  possible_mothers <- table(parents[2:nrow(parents),i])
  print("Här är i")
  print(i)
  print(possible_mothers)
  if(length(possible_mothers)>1){
    print("Long")
    x <- possible_mothers[names(possible_mothers)!="0"]
    x <- x[names(x)!="NP"]
    print(x)
    print(names(x))
    if(max(x)==1){
      parents_new[1,i] <- min(names(x))
    }else{
      parents_new[1,i] <- names(x)[1]
    }
  }else{
    possible_mothers <- table(parents[2:nrow(parents),i])
    parents_new[1,i] <- names(possible_mothers)[1]
  }
  i <- i+1
}

parents_new <- as.vector(as.numeric(parents_new))

frac.table_new[,4] <- as.matrix(as.numeric(frac.table_new[,4])*0.1)
frac.table_new[,7] <- as.matrix(as.numeric(frac.table_new[,7])*0.1)

fish = createFishObject(frac.table_new,parents_new,timepoints=timepoints,fix.missing.clones = TRUE)
fish = layoutClones(fish)

library(viridis)
library(scales)
#show_col(viridis_pal()(17))
fish <- setCol(fish,viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1,option = "D")(19))
fishPlot(fish,shape="spline",title.btm="LUNOV",
         cex.title=0.5, vlines=timepoints,
         vlab=c("Day 0","Day 100","Day 150","","Day 250","Day 300","","Day 400","Day 600"),bg.type="none")
```
