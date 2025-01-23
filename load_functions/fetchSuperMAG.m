%{
fetchSuperMAG.m
SuperMAG inventory, data and indices fetching script
Author: S. Antunes
May 2021, v 1.0

(c) 2021  The Johns Hopkins University Applied Physics Laboratory LLC.
All Rights Reserved. 

This material may be only be used, modified, or reproduced by or for
the U.S. Government pursuant to the license rights granted under the
clauses at DFARS 252.227-7013/7014 or FAR 52.227-14. For any other
permission, please contact the Office of Technology Transfer at
JHU/APL.

NO WARRANTY, NO LIABILITY. THIS MATERIAL IS PROVIDED "AS IS." JHU/
APL MAKES NO REPRESENTATION OR WARRANTY WITH RESPECT TO THE
PERFORMANCE OF THE MATERIALS, INCLUDING THEIR SAFETY, EFFECTIVENESS,
OR COMMERCIAL VIABILITY, AND DISCLAIMS ALL WARRANTIES IN THE
MATERIAL, WHETHER EXPRESS OR IMPLIED, INCLUDING (BUT NOT LIMITED TO)
ANY AND ALL IMPLIED WARRANTIES OF PERFORMANCE, MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF
INTELLECTUAL PROPERTY OR OTHER THIRD PARTY RIGHTS. ANY USER OF THE
MATERIAL ASSUMES THE ENTIRERISK AND LIABILITY FOR USING THE
MATERIAL. IN NO EVENT SHALL JHU/APL BE LIABLE TO ANY USER OF THE
MATERIAL FOR ANY ACTUAL, INDIRECT, CONSEQUENTIAL, SPECIAL OR OTHER
DAMAGES ARISING FROM THE USE OF, OR INABILITY TO USE, THE MATERIAL,
INCLUDING, BUT NOT LIMITED TO, ANY DAMAGES FOR LOST PROFITS.
%}

% Basic usage:
% for inventory: mydata=fetchSuperMAG('inventory',userid,start,extent)
% for indices:   mydata=fetchSuperMAG('indices',userid,start,extent,flags)
% for data:      mydata=fetchSuperMAG('data',userid,start,extent,flags,station)
%
% Test cases: assumes userid set to your SuperMAG access id
% mydata=fetchSuperMAG('inventory',userid, '2019-11-15T10:40', 3600);  
% mydata=fetchSuperMAG('indices',userid, '2019-11-15T10:40', 3600,'plusall,noisy');
% mydata=fetchSuperMAG('data',userid, '2019-11-15T10:40', 3600,'sunall,noisy');
% ========================================================================

function SuperMAG_data = fetchSuperMAG(category,userid,start,extent,flags,station)
  % do some initial checking of input variables (and init the data return)
  SuperMAG_data='';
  % choices are to fetch 'data' or 'indices' or 'inventory'
  category=lower(category);
  if (strcmp(category,'data') | strcmp(category,'indices') | strcmp(category,'inventory') ) == 0
    error="Error: Please state whether you want 'inventory', 'data' or 'indices'";
    disp(error)
    SuperMAG_data=error;
    return
  end

  if nargin < 4
    error="Please provide category, userid, start and extent";
    disp(error)
    SuperMAG_data=error;
    return
  end
  
  if strcmp(category,'data') & nargin < 6
    error="Please include a station for fetching data";
    disp(error)
    SuperMAG_data=error;
    return
  end

  % DETAILS
  % For category='data', always returns:
  %        [list of stations]
  %

  
  % For category='data', always returns:
  %        tval, n, e, z
  % additional flags to return data include:
  %        mlt,mag,geo,decl,sza
  %        (or 'all')
  % e.g. either 'mlt,mag' or ["mlt mag"] will work
  % optional flags that modify data are:
  %        delta=start
  %        baseline=none or baseline=yearly
  %

  % For category='indices', always returns:
  %        tval
  % additional flags to return data include:
  %        indicesall (or its alias: all)
  %  (or any of)
  %        baseall, sunall, darkall, regionalall, plusall
  %  (or specify individual items to include, from the sets below)
  %        
  %basekeys=["sme","sml","smu","mlat","mlt","glat","glon","stid","num"];
    % sunkeys: alias allowed of SUN___ -> ___s
  %sunkeys=["smes","smls","smus","mlats","mlts","glats","glons","stids","nums"];
    % darkkeys: alias allowed of DARK___ -> ___d
  %darkkeys=["smed","smld","smud","mlatd","mltd","glatd","glond","stidd","num"];
    % regkeys: alias allowed of REGIONAL___ -> ___r
  %regkeys=["smer","smlr","smur","mlatr","mltr","glatr","glonr","stidr","numr"];
  %pluskeys=["smr","ltsmr","ltnum","nsmr"];
    
% ========================================================================
  % Notes on this implementation:
  % This code has 3 sections:
  % PART A: Here is the ingest of user-provide arguments
  % PART B: Here is the actual server code
  % PART C: Here is optional usage of the data (flag 'noisy' will run this)
  %
  % testurlD='https://supermag.jhuapl.edu/services/data-api.php?fmt=json&nohead&logon=YOURNAME&start=2019-10-15T10:40&extent=3600&mlt&mag&geo&decl&sza' 
  % testurlI='/supermag.jhuapl.edu/services/indices.php?fmt=json&logon=YOURNAME&start=2019-11-02T20:24&extent=3600&indices=all'
  % testurlI='/supermag.jhuapl.edu/services/inventory.php?fmt=json&logon=YOURNAME&start=2019-11-02T20:24&extent=3600&station=HBK'
  
  % ('json' + 'nohead' tells SuperMAG server to send a proper JSON objects)
  % ('matlab' flag is a courtesy to help SuperMAG team track stats on usage)
  %
  % note Matlab has an 'arguments' block for parsing on for R2019b or later
  % so we do our own intelligent argument parsing.  Works for Matlab >2017

% ========================================================================
  % PART A: Here is the ingest of user-provide arguments
  %
  
  noisy=0; % initially stay quiet, set flags of 'noise' if you want diagnostics

  % accepts start either as a string (unparsed), a datatime, or an array
  % of format [YYYY MO DD HH MM SS] (SS is optional)
  % e.g. '2019-10-15T10:40', "2019-10-15T10:40", 2019-10-15T10:40:00 datetime,
  % or [2019 10 15 10 40 0] are all valid inputs

  if numel(start) == 5 | numel(start) == 6 | isdatetime(start)
    % user chose to give us a datetime object or array, so convert it
    disp("formatting")
    mystart = string(datetime(start,'Format','yyyy-MM-dd''T''HH:mm:00'))
  else
    % hoping they gave us the proper formatted datetime string
    mystart=start;
  end

  mystation='';

  if strcmp(category,'inventory')
    page='inventory.php';
  elseif strcmp(category,'data')
    page='data-api.php';
    % only 'data' needs a specific station, 'indices' is gridded instead
    mystation=append('&station=',station);
    % full list of additional data sets available
    datakeys=["mlt","mag","geo","decl","sza","baseline=none","baseline=yearly","delta=start"];
  elseif strcmp(category,'indices')
    page='indices.php';
    % note 'indices=all' is datakeys + sunkeys + darkkeys + regkeys + pluskeys
    basekeys=["sme","sml","smu","mlat","mlt","glat","glon","stid","num"];
    % sunkeys: alias allowed of SUN___ -> ___s
    sunkeys=["smes","smls","smus","mlats","mlts","glats","glons","stids","nums"];
    % darkkeys: alias allowed of DARK___ -> ___d
    darkkeys=["smed","smld","smud","mlatd","mltd","glatd","glond","stidd","numd"];
    % regkeys: alias allowed of REGIONAL___ -> ___r
    regkeys=["smer","smlr","smur","mlatr","mltr","glatr","glonr","stidr","numr"];
    pluskeys=["smr","smrlt","ltnum","nsmr"];
    indiceskeys=[basekeys sunkeys darkkeys regkeys pluskeys];
    % 'all' means all the above
    
    imfkeys=["bgse","bgsm","vgse","vgsm"]; % or imfall for all these
    swikeys=["pdyn","epsilon","newell","clockgse","clockgsm","density"]; % or swiall for all these

  end

  myflags='';
  indices='&indices=';
  swi='&swi=';
  imf='&imf=';

  if nargin > 4
    % for 'data' and 'indices' (or 'inventory' if 'noisy' flag added)--

    % parse flags, any of the following
    % flags can either be a string (e.g. 'mlt,mag,baseline=yearly')
    % or a list (e.g.a ["mlt" "mag" "baseline=yearly"]

    % pre-parse, if string, into a list so we can walk through it
    if ischar(flags)
      flags=strtrim(string(strsplit(flags,',')));
    end

    % match any 'known' flags the user requested and add to url string
    for i=1:numel(flags)
      chk=lower(strtrim(flags(i)));

      % special keys
      if strcmp(chk,'noisy'); noisy=1; end;

      if strcmp(category,'data')
        % check for 'all', also individual keys, and assemble url flags
        if strcmp(chk,'all')
	  myflags=append(myflags,'&mlt&mag&geo&decl&sza');
	end
	for ikey=1:length(datakeys)
	  if strcmp(chk,datakeys(ikey))
	    myflags=append(myflags,'&',datakeys(ikey));
	  end
	end
      end

      if strcmp(category,'indices')

        % check for the '*all', also individual keys, and assemble url flags
	if strcmp(chk,'all'); indices=append(indices,'all,'); end;
	if strcmp(chk,'indicesall'); indices=append(indices,'all,'); end;
	if strcmp(chk,'imfall'); imf=append(imf,'all,'); end;
	if strcmp(chk,'swiall'); swi=append(swi,'all,'); end;
        % available keywords, we allow both the url version and the
        % aliases of "SUN___ -> ___s", "DARK___ -> ___d", "REGIONAL___ -> ___r"
	for ikey=1:length(indiceskeys)
	  mykey=indiceskeys(ikey);        % the legit url keyword to check
	  sunkey=append("sun",mykey);       % allow for alias
	  darkkey=append("dark",mykey);     % allow for alias
	  regkey1=append("regional",mykey);  % allow for alias
	  regkey2=append("reg",mykey);  % allow for alias
	  if strcmp(chk,mykey)
	    indices=append(indices,mykey,','); % base key is correct
	  elseif strcmp(chk,sunkey)
	    indices=append(indices,mykey,'s,'); % alias, so base key + 's'
	  elseif strcmp(chk,darkkey)
	    indices=append(indices,mykey,'d,'); % alias, so base key + 'd'
	  elseif (strcmp(chk,regkey1) | strcmp(chk,regkey2))
	    indices=append(indices,mykey,'r,'); % alias, so base key + 'r'
	  end
	end

	for ikey=1:length(swikeys)
	  if strcmp(chk,swikeys(ikey)); swi=append(swi,swikeys(ikey),','); end;
	end
	for ikey=1:length(imfkeys)
	  if strcmp(chk,imfkeys(ikey)); imf=append(imf,imfkeys(ikey),','); end;
	end

        % more aliases for the user
	if strcmp(chk,'baseall'); indices=append(indices,join(basekeys,','),','); end;
	if strcmp(chk,'sunall'); indices=append(indices,join(sunkeys,','),','); end;
	if strcmp(chk,'darkall'); indices=append(indices,join(darkkeys,','),','); end;
	if (strcmp(chk,'regionalall') | strcmp(chk,'regall')); indices=append(indices,join(regkeys,','),','); end;
	if strcmp(chk,'plusall'); indices=append(indices,join(pluskeys,','),','); end;
	
      end

    end
  end

  if strcmp(category,'indices')
    % clean it up a bit by removing extraneous tags/characters
    if strcmp(indices,"&indices="); indices=""; end;
    if strcmp(swi,"&swi="); swi=""; end;
    if strcmp(imf,"&imf="); imf=""; end;
    % add them together
    myflags=append(indices,swi,imf);
    % a little more cleaning for tidiness, removes extraneous commas
    myflags=strrep(myflags,',&','&');
    myflags=regexprep(myflags,',$','');
  end
  
% ========================================================================
  % PART B: Here is the actual server code
  % now construct the full URL for the data fetch  
  SERVER='https://supermag.jhuapl.edu/services/';

  logon=append('&logon=',userid);

  spec='?matlab&nohead';
  mystart=append('&start=',mystart);
  extent=append('&extent=',string(extent));

  % inventory uses SERVER/page/logon/start/extent
  % data-api uses above + station/flags
  % indices ditto sans station
  url = append(SERVER,page,spec,logon,mystart,extent,myflags,mystation);

  % optional, if 'noisy' provided then displays url
  if noisy; disp(url); end; 
  
  if strcmp(category,'inventory')
    opts = weboptions('ContentType','text','Timeout',Inf);
  else
    opts = weboptions('ContentType','json','Timeout',Inf);
  end
  
  try
    SuperMAG_data = webread(url,opts);
  catch err
    error=err.message;
    disp(error)
    SuperMAG_data=error;
    return
  end

  if strcmp(category,'inventory')
    % 'inventory' is not JSON data so we array-ize it
    SuperMAG_data = strsplit(SuperMAG_data,'\n');
    SuperMAG_data=SuperMAG_data(2:end-1) % leading field was # of stations
  end

  % Done!  The rest of this code is extra, including examples of using the data
% ========================================================================
  % PART C: Here is optional usage of the data (flag 'noisy' will run this)
  if noisy & strcmp(category,'inventory')
    for i = 1:length(SuperMAG_data)
      disp(SuperMAG_data{i})
    end
  elseif noisy
    % data is a structure, isstruct(data)
    % optional, if 'noisy' provided then echos a subset of rows and their
    % values to screen as well
    % Also useful as example code on how to walk the SuperMAG data
    % note that arrays start at 1, not 0
    if noisy
      disp("Printing 3 sample rows of data, for validation")
      for irow=1:min(3,length(SuperMAG_data))
	fprintf("\nRow %d\n",irow)
	keys=fieldnames(SuperMAG_data);
	for i = 1:length(keys)
	  key=keys{i};
	  value=getfield(SuperMAG_data,key);
	  if isstruct(value)
	    % item is itself a structure of data items
	    subkeys=fieldnames(value);
	    for subi = 1:length(subkeys)
	      subkey=subkeys{subi};
	      subvalue=getfield(value,subkey);
	      disp(key+"."+subkey+":"+subvalue);
	    end
	  else
	    % item is a single value (float or string)
	    disp(key+":"+value);
	  end
	end
        % data(1).N = itself a struct with fields nez, geo, so...
        % disp(mycell.N.geo) etc works
      end
    end

  end

%{
sample fetches to indicate the variable names returned:

========================================================================
  sample fetch of data
  data = 
    60×1 struct array with fields:
      tval
      ext
      iaga
      glon
      glat
      mlon
      mlat
      mlt
      mcolat
      decl
      sza
      N
      E
      Z

   e.g. data(1) gives 1 cell, data(1).tval gives that data item
   length(data) is 60, so here is one sample plotting of 2 data items:

sm_data=fetchSuperMAG('data',userid,'2019-11-15T10:40',3600,'all,delta=start,baseline=yearly','HBK');

% simple plot of two nested structure elements
tval=[sm_data.tval]
% extract the N.geo and N.nez elements
N=[sm_data.N]   
% plot the N.geo and N.nez elements
plot( [tval], [N.geo])
hold on
plot( [tval], [N.nez])
hold off
title("N geo vs N nez")
xlabel("time")
ylabel("N vector")
%
% Alternate way to extract a nested structure element such as sm_data.N.geo
N_geo = arrayfun(@(k) sm_data(k).N.geo, 1:numel(sm_data));
plot([sm_data.tval], N_geo)
%
% if you prefer everything as arrays instead of structures:
TVAL=sm_data[*].tval
N_NEZ = arrayfun(@(k) sm_data(k).N.nez, 1:numel(sm_data));
E_NEZ = arrayfun(@(k) sm_data(k).E.nez, 1:numel(sm_data));
Z_NEZ = arrayfun(@(k) sm_data(k).Z.nez, 1:numel(sm_data));
MLT=[sm_data.mlt];
MCOLAT=[sm_data.mcolat];
MLON=[sm_data.mlon];
MLAT=[sm_data.mlat];
GLON=[sm_data.glon];
GLAT=[sm_data.glat];
SZA=[sm_data.sza];

========================================================================
   sample fetch of inventory:
   data =
    1x174 cell array


   sample fetch of indices:
   data = 
  60x1 struct array with fields:
    tval
    SME
    SML
    SMLmlat
    SMLmlt
    SMLglat
    SMLglon
    SMLstid
    SMU
    SMUmlat
    SMUmlt
    SMUglat
    SMUglon
    SMUstid
    SMEnum
    SMEs
    SMLs
    SMLsmlat
    SMLsmlt
    SMLsglat
    SMLsglon
    SMLsstid
    SMUs
    SMUsmlat
    SMUsmlt
    SMUsglat
    SMUsglon
    SMUsstid
    SMEsnum
    SMEd
    SMLd
    SMLdmlat
    SMLdmlt
    SMLdglat
    SMLdglon
    SMLdstid
    SMUd
    SMUdmlat
    SMUdmlt
    SMUdglat
    SMUdglon
    SMUdstid
    SMEdnum
    SMEr
    SMLr
    SMLrmlat
    SMLrmlt
    SMLrglat
    SMLrglon
    SMLrstid
    SMUr
    SMUrmlat
    SMUrmlt
    SMUrglat
    SMUrglon
    SMUrstid
    SMErnum
    smr
    smr00
    smr06
    smr12
    smr18
    smrnum
    smrnum00
    smrnum06
    smrnum12
    smrnum18

sm_data=fetchSuperMAG('indices',userid,'2019-11-15T10:40',3600,'all');

; plot time vs a data item
plot([sm_data.tval],[sm_data.SMEs])
hold on
plot([sm_data.tval],[sm_data.SMLs])
hold off


% plot time vs a multi-dimensional data item
tval=[sm_data.tval];
y=[sm_data.SMLr];
nrows =	length(tval);
hours = 0:23;
for i=1 : nrows
  plot( hours, y(:,i) )
  hold on
end
hold off
title('Sample Plot, SMLr vs Hour across multiple days') ;



========================================================================
%}

end
