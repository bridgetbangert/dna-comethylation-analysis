args <- commandArgs(TRUE)
setwd("/home/bmb191/math5376D/Final.Project/Output.Files/Correlation.Counts")

# investigation to get the difference between 2 correlation count arrays

# load in data
cnt_array_1 <- read.table(args[1], header=F)
cnt_array_2 <- read.table(args[2], header=F)

# function to create df comparing the two, write it to text file
create_diff_df <- function(arr1, arr2) {
	count_df <- data.frame(
		Counts_1=arr1,
		Counts_2=arr2
	)
	cbind(1:nrow(count_df), count_df)
}

# write diff df to leap server
write_diff_df <- function(arr1, arr2) {
  	write.table(
  		create_diff_df(arr1, arr2), 
  		args[3], 
  		row.names=FALSE, 
  		col.names=FALSE
	)
}

# write diff df to leap server
find_diff_cnts <- function(arr1, arr2) {
	count_df <- create_diff_df(arr1, arr2)
  	count_df[count_df$Counts_1 != count_df$Counts_2, ]
}

write_diff_df(cnt_array_1, cnt_array_2)

