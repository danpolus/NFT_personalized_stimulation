% runs NFTsim from every directory
% inputs:
%   configuration - file name
%   nftsim_path -
%   cd_flg - move to nft folder or stay in the current one
% outputs:
%   nf_struct
%
function nf_struct = run_nftsim_pwd(configuration, nftsim_path, cd_flg)

if cd_flg
    fp = '.';
    addpath(nftsim_path);
    work_dir = pwd;
    cd(nftsim_path);
else
    fp = nftsim_path;
end

in_fn = [fp '\configs\' configuration '.conf'];
out_fn = [fp '\output\' configuration '.output'];
disp(['NFTsim data generation: ' configuration]);
[status, cmdout] = system([fp '\bin\nftsim.exe -i ' in_fn ' -o ' out_fn]);
if status == 0
    disp(cmdout);
    nf_struct = nf.read(out_fn);
    nf_struct.conf_file = configuration;
else
    error(cmdout);
end

if cd_flg
    cd(work_dir);
end
