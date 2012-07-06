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
				from trackdb""")]

def read_tracks_available(dbpath):
	return [r for (r,) in 
		execute(dbpath, """select name from sqlite_master
				where type = 'table' and name not in 
				('chromInfo', 'trackdb')""")]

queries = {
	"bed 3": """select chrom, chromStart, chromEnd from {}""",
	"bed 6": """select chrom, chromStart, chromEnd, 
		name, score, strand from {}"""}


def write_bed(dbpath, track, outpath):
	with closing(sqlite3.connect(dbpath)) as conn:
		with closing(conn.cursor()) as c:
			#TODO: add more queries
			if not track.type.startswith("bed"):
				return
			try:
				type = "bed 6" if int(track.type[4]) > 6 else "bed 3"
			except ValueError:
				type = "bed 3"
			qry = queries.get(type, queries["bed 3"])
			print track.name, track.type, type
			c.execute(qry.format(track.name))
			outdir = os.path.dirname(outpath)
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			with open(outpath, "w") as bed:
				for row in c:
					bed.write("\t".join(map(str,row))+"\n")

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
					write_bed(dbpath, trk, path)


if __name__ == "__main__":
	if len(sys.argv) != 2:
		raise SystemExit("USAGE: python -m genome_runner.unpack_db")
	dump_trackdb(sys.argv[1], "data")
