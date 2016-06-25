import os
import logging
from logging import FileHandler,StreamHandler
import subprocess
import collections

logger = logging.getLogger('genomerunner.server')

def retrieve_files(up_files, outputdir, id):
	"""
	Used to retrieve either a single file or multiple files uploaded by the user
	:param up_files: the files that are returned by the user
	:param outputdir: The directory to which the file should be uploaded to sever
	:param file_list_path: File into which the full path of file on server is written (i.e. '/path/.fois')
	:param id: the id of the grsnp job
	:return: List of all of paths of the uploaded files
	"""
	uploaded_files = []
	if not isinstance(up_files, (list)):
		up_files = [up_files]  # makes a list if only one file uploaded
	if up_files[0] and up_files[0].filename != "":
		for f in up_files:
			f_filename = f.filename
			f_outpath = os.path.join(outputdir,f_filename)
			extension = f_filename.split('.')[-1]
			if not os.path.exists(f_outpath):
				with open(f_outpath,'wb') as out:
					if f != None and f.filename != "":
						logger.info("Received uploaded file (name={}, id={})".format(f.filename,id))
						while True:
							data = f.file.read(8192)
							# TODO find empty lines
							# data = os.linesep.join([s for s in data.splitlines() if s ])

							if not data:
								break
							out.write(data)


				if extension not in ["gz", "bb"]:
					# use dos2unix to remove \r from end of lines
					out = subprocess.Popen(['dos2unix {}'.format(f_outpath)], shell=True, stdout=subprocess.PIPE,
										   stderr=subprocess.PIPE)
					out.wait()
			else:
				logger.error("id={} Upload file already exists at {}".format(id, f_outpath))
				print "id={} Upload file already exists at {}".format(id, f_outpath)
			uploaded_files.append((f_outpath))
	return uploaded_files

def retrieve_text(up_text, out_file_path, id):
	"""
	Used to retrieve the bed data uploaded via text box.
	:param up_text: the bed text retured by the user
	:param outputdir: The file to which the bed text should be uploaded to sever
	:param file_list_path: File into which the full path of file on server is written (i.e. '/path/.fois')
	:param id: the id of the grsnp job
	:return: path to the bed file on the server in list form [path]
	"""
	uploaded_t_p = []
	try:
		if up_text != None and up_text != "":
			with open(out_file_path, "wb") as out:
				logger.info('Received raw text data (name = {}, id={})'.format(os.path.split(out_file_path)[-1], id))
				data = os.linesep.join([s for s in up_text.split("\n") if s != ""])
				out.write(data)
			# use dos2unix to remove \r from end of lines
			out = subprocess.Popen(['dos2unix {}'.format(out_file_path)], shell=True, stdout=subprocess.PIPE,
								   stderr=subprocess.PIPE)
			out.wait()
			uploaded_t_p = [out_file_path]

	except Exception, e:
		logger.error("id={}".format(id) + str(e))
		return "ERROR: Unable to process custom bed file"
	return uploaded_t_p

def retrieve_group(up_group_folder):
	ls_grp = []
	if up_group_folder and up_group_folder != "":
		ls_grp = [os.path.join(up_group_folder, f) for f in os.listdir(up_group_folder) if
				  os.path.isfile(os.path.join(up_group_folder, f)) and not f.endswith(".tbi")]
	return ls_grp


def purge_duplicate_features(file_path_list):
	'''
	Removes all filespaths from the buttom up that have the same basename.
	[dir/file1.txt, dir/file2.bed, dir/file1.bed, dir/file2.bed.gz] ==> [dir/file1.txt, dir/file2.bed]
	:param file_path_list: list of file paths to filter
	:return:
	'''
	list_base =  [base_name(x) for x in file_path_list]
	dups = collections.defaultdict(list)
	for i, e in enumerate(list_base):
		dups[e].append(i)
	indices_delete = []
	for k, v in sorted(dups.iteritems()):
		if len(v) >= 2:
			indices_delete = indices_delete + v[1:]

	# delete all duplicate indeces in reverse order
	return [i for j, i in enumerate(file_path_list) if j not in indices_delete]

def base_name(k):
	return os.path.basename(k).split(".")[0]
