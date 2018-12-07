# convert factor to numeric
fac2num <-
    function(v)
    {
        nam <- names(v)
        v <- as.numeric(as.character(v))
        names(v) <- nam
        v
    }

# example
x <- factor(c(5, 4, 1, 4, 5))
as.numeric(x)  # doesn't do what we want
fac2num(x)
