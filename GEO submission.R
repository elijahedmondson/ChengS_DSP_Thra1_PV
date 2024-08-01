

file_name <- "C:/Users/edmondsonef/Desktop/230816_VH00687_252_AAC5TGFHV.tar"
untar(file_name)
getwd()






DCCFiles <- list.files("C:/Users/edmondsonef/Desktop/GEO/processed_data/", full.names=TRUE)

## compute check sum for every file in package
x <- tools::md5sum(DCCFiles)

## exclude any existing "MD5" file
x <- x[names(x) != "MD5"]
outDir <- "C:/Users/edmondsonef/Desktop/"
## write each result out to the "MD5" file
cat(paste(x, names(x), sep = " *"), sep = "\n", 
    file = file.path(outDir, "MD5"))


write.csv(x, "C:/Users/edmondsonef/Desktop/data.csv")
