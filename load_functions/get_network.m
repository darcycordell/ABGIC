function [L,S,T] = get_network
%%
get_files;

curdir = pwd;

cd(networkFilePath)
load(lineFile)
load(subFile)
load(tranFile)
cd(curdir)

if ~isfield(L,'fromBus')
    %If there are no buses specified in LINES.mat file, then we compute
    %them using some defaults:
    %       (1) Assume a single transformer in each substation
    %       (2) Assume 240 kV transformers to terminations are wye-delta
    %       (3) Assume 500 kV to 240 kV transformers are autotransformers
    [L,T] = get_buses(L,S);
end

segLength = 5000; %downsample lines so that they are divided into x meters
L = get_downsampled_line(L,segLength);


%Get line lengths. Note that this script also flips line coordinates so
%that the line coordinates begin at the "from bus" and go to the "to bus"
L = calc_line_length(L,S);

if ~isfield(L,'Resistance')
    %If resistances have not been specified in the LINES.mat file, then
    %compute them using Horton defaults
    L = calc_line_resistance(L);
end