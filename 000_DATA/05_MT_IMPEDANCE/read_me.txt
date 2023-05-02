MT DATA NOTES:

AB_BC_MT_DATA_512_sites.mat:
	This is the original data file given to me by Cedar in 2021. This is the file that was
	used for the 2021 Space Weather paper. This data was edited by Cedar and Q/C'd by me.
	It was also the dataset used for inversion modelling.

AB_BC_MT_DATA_527_sites_230310.mat:
	This is the same as the above but includes an additional 15 sites collected by Unsworth
	group in summer of 2022. These sites begin with prefix "AER" and are mostly collected
	in the area west of Fort McMurray and east of Red Earth Creek. By coincidence, these
	data fill an important gap for geoelectric field modelling since this is where the new
	500 kV transmission line from Fort McMurray to Genesee is located.

AB_BC_MT_DATA_526_sites_230321.mat:
	This is the same as AB_BC_MT_DATA_527 except with sab235 removed. SAB235 is a site located
	in eastern Alberta near Provost on the SK border. This site has always appeared as a bit of
	an outlier with larger geoelectric fields than adjacent sites. I never bothered to edit it
	since it is not near any transmission lines and thus doesn't affect the line voltage 
	estimates. However, when looking at peak geoelectric fields for the whole region