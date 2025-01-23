%Download data for selected sites directly from INTERMAGNET using URL call

%see here: https://imag-data.bgs.ac.uk/GIN_V1/WebServiceURLGenerator.html

sites = {'NEW','BOU','YKC','MEA','FCC','VIC','BRD','SIT'};
dates = {'2024-05-10','2024-05-11'};

for i = 1:length(sites)
    for j = 1:length(dates)

    url = ['https://imag-data.bgs.ac.uk/GIN_V1/GINServices?Request=GetData&observatoryIagaCode=',sites{i},'&samplesPerDay=second&dataStartDate=',dates{j},'&dataDuration=1&publicationState=reported&orientation=XYZF&format=iaga2002'];

    data = webread(url);

    filename_out =   [lower(sites{i}),erase(dates{j},'-'),'vsec.sec'];      
    fid = fopen(filename_out, 'w' ); %// open file to writing
    fprintf( fid, data ); %// print string to file
    fclose( fid ); %// don't forget to close the file

    end
    sites{i}
end

