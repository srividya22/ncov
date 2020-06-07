#!/usr/bin/env bash

# Usage: drophdrs [data [exclude]]

awk 'BEGIN{FS = OFS = ","}
dbg{    printf("FILENAME=%s, FNR=%d, NR=%d, NF=%d, $0=\"%s\"\n",
                FILENAME, FNR, NR, NF, $0)
}
FNR==1{ if(dbg) printf("%s file header with %d fields: %s\n", FILENAME, NF, $0)
        if(FNR==NR) {
                efn = FILENAME # Save filename of exclude file for diagnostics.
                next
        }
        # Determine which fields to skip from headers in the data file.
        for(i = 1; i <= NF; i++) if($i in skiphdr) {
                sf[i]
                delete skiphdr[$i]
                if(dbg) printf("Field %d added to sf[] for header %s.\n", i, $i)
        }
        first = 1
        for(i in skiphdr) {
                if(first) {
                        first = 0
                        printf("File %s will not be processed because:\n",
                                FILENAME)
                }
                printf("\theader \"%s\" in exclude file (%s) was not found\n",
                        i, efn, FILENAME)
        }
        if(first == 0) exit 1
}
FNR==NR{# gather names of columns to be skipped from exclude (1st) file
        skiphdr[$1]
        if(dbg) printf("%s added to skiphdr\n", $1)
        next
}
{       sep = ""
        for(i = 1; i <= NF; i++)
                if(!(i in sf)) {
                        printf("%s%s", sep, $i)
                        sep = OFS
                }
        printf("\n")
}' ${2:-exclude} ${1:-input}
