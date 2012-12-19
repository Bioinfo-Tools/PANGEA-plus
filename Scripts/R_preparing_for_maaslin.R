table10f <- read.table("10f_functions4.txt", header=F)
table10f <- read.table("10f_functions4.txt", header=F)
table11b <- read.table("11b_functions4.txt", header=F)
table18j <- read.table("18j_functions4.txt", header=F)
table20g <- read.table("20g_functions4.txt", header=F)
table24d <- read.table("24d_functions4.txt", header=F)
table27e <- read.table("27e_functions4.txt", header=F)
table29f <- read.table("29f_functions4.txt", header=F)
table32f <- read.table("32f_functions4.txt", header=F)
table34n <- read.table("34n_functions4.txt", header=F)
table3n <- read.table("3n_functions4.txt", header=F)
table44j <- read.table("44j_functions4.txt", header=F)
table4j <- read.table("4j_functions4.txt", header=F)
table6f <- read.table("6f_functions4.txt", header=F)
table9l <- read.table("9l_functions4.txt", header=F)


require(data.table)

require(zoo)

frame10f <- data.frame(table10f)
frame11b <- data.frame(table11b)
frame18j <- data.frame(table18j)
frame20g <- data.frame(table20g)
frame24d <- data.frame(table24d)
frame27e <- data.frame(table27e)
frame29f <- data.frame(table29f)
frame32f <- data.frame(table32f)
frame34n <- data.frame(table34n)
frame3n <- data.frame(table3n)
frame44j <- data.frame(table44j)
frame4j <- data.frame(table4j)
frame6f <- data.frame(table6f)
frame9l <- data.frame(table9l)

test1 <- merge.data.frame(frame10f, frame11b, all=TRUE, fill=c("0"), by=c("V1"))
test2 <- merge.data.frame(test1, frame18j, all=TRUE, fill=c("0"), by=c("V1"))
test3 <- merge.data.frame(test2, frame20g, all.x=TRUE, fill=c("0"), by=c("V1"))
test4 <- merge.data.frame(test3, frame24d, all.x=TRUE, fill=c("0"), by=c("V1"))
test5 <- merge.data.frame(test4, frame27e, all.x=TRUE, fill=c("0"), by=c("V1"))
test6 <- merge.data.frame(test5, frame29f, all.x=TRUE, fill=c("0"), by=c("V1"))
test7 <- merge.data.frame(test6, frame32f, all.x=TRUE, fill=c("0"), by=c("V1"))
test8 <- merge.data.frame(test7, frame34n, all.x=TRUE, fill=c("0"), by=c("V1"))
test9 <- merge.data.frame(test8, frame3n, all.x=TRUE, fill=c("0"), by=c("V1"))
test10 <- merge.data.frame(test9, frame44j, all.x=TRUE, fill=c("0"), by=c("V1"))
test11 <- merge.data.frame(test10, frame4j, all.x=TRUE, fill=c("0"), by=c("V1"))
test12 <- merge.data.frame(test11, frame6f, all.x=TRUE, fill=c("0"), by=c("V1"))
test13 <- merge.data.frame(test12, frame9l, all.x=TRUE, fill=c("0"), by=c("V1"))

sink("outfile.txt")
print(test13)
write.table(test13, file = "outfile.txt", append = FALSE, quote = FALSE, sep = "\t",
eol = "\n", na = "0", dec = ".", 
col.names = TRUE, qmethod = c("escape", "double"),
row.names = FALSE,
fileEncoding = "")

sink()