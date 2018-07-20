#!/usr/bin/awk -f
{
    sqFile=sample"_sq.txt"
    pgFile=sample"_pg.txt"
    if ($1 ~ /^@SQ/) { print $0 > sqFile }
    if ($1 ~ /^@PG/) { print $0" sample: "fileBase > pgFile }
}
