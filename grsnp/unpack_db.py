#!/usr/bin/env pypy
import sys, sqlite3, os
from collections import namedtuple
from logging import debug
from contextlib import closing
from itertools import groupby
from operator import attrgetter

Track = namedtuple("Track", 
		["name", "type", "grp", "visibility", "priority", "html"])

def execute(dbpath, qry):
	with closing(sqlite3.connect(dbpath)) as conn:
		with closing(conn.cursor()) as c:
			c.execute(qry)
			return [r for r in c]
		
def rsgroupby(k, iterable):
	return groupby(sorted(iterable, key=k, reverse=True), k)

def read_trackdb(dbpath):
	return [Track(*r) for r in execute(dbpath,
		"""select tableName,type,grp,
			visibility,priority,html 
				from trackDb""")]

def read_tracks_available(dbpath):
	return [r for (r,) in 
		execute(dbpath, """select name from sqlite_master
				where type = 'table' and name not in 
				('chromInfo', 'trackDb')""")]

def write_rows(path, rows):
    outdir = os.path.dirname(outpath)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    with open(outpath, "w") as bed:
        for row in c:
            bed.write("\t".join(map(str,row))+"\n")


def write_bed(dbpath, track, outpath):
    try:
        type = "bed 6" if int(track.type[4]) > 6 else "bed 3"
    except ValueError:
        type = "bed 3"
   
    if type == "bed 3":
        qry = """select chrom, chromStart, chromEnd from {}"""
    else:
        qry = """select chrom, chromStart, chromEnd, 
		name, score, strand from {}"""
    write_rows(outpath, execute(dbpath, qry.format(track.name)))

def write_gene_bed(dbpath,track,outpath):
    rows = execute(dbpath,"""select chrom,txStart,txEnd,exonStarts,exonEnds,name from {}""".format(track.name))
    dir = os.path.dirname(outpath)
    genes = [(chrom, txStart, txEnd, name) for (chrom,txStart,txEnd,name,_,_) in rows]
    write_rows(genes, os.path.join(dir,track.name+".bed"))
    exons = []
    for (chrom,_,_,exonStarts,exonEnds,name) in rows:
        for (s,e) in zip(exonStarts.split(","), exonEnds.split(",")):
            exons.append((chrom,s,e,name))
    write_rows(exons, os.path.join(dir,track.name+"_exons.bed"))

def handle_track(dbpath, track, outpath):
    print "writing file"
    if track.type.startswith("bed"):
        write_bed(dbpath, track, outpath)
    elif track.type.startswith("genePred"):
        write_gene_bed(dbpath,track,outpath)

def dump_trackdb(dbpath, outdir):
	available = set(read_tracks_available(dbpath))
	trackdb = read_trackdb(dbpath)
	for grp, grp_trks in rsgroupby(attrgetter("grp"), trackdb):
		for vis, vis_trks in rsgroupby(attrgetter("visibility"), grp_trks):
			vis = "Tier"+str(vis)
			for trk in vis_trks:
				if trk.name in available:
					path = os.path.join(outdir,grp,vis,trk.name+".bed")
					debug(trk.name)
					handle_track(dbpath, trk, path)


if __name__ == "__main__":
	if len(sys.argv) != 2:
		raise SystemExit("USAGE: python -m genome_runner.unpack_db")
	dump_trackdb(sys.argv[1], "data")
