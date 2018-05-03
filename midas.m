%  MIDAS
%  Version 1.0.0 --- January 2018
%
%  Section of Biomedical Image Analysis
%  Department of Radiology
%  University of Pennsylvania
%  Richard Building
%  3700 Hamilton Walk, 7th Floor
%  Philadelphia, PA 19104
%
%  Web:   https://www.med.upenn.edu/sbia/
%  Email: sbia-software at uphs.upenn.edu
%
%  Copyright (c) 2018 University of Pennsylvania. All rights reserved.
%  See https://www.med.upenn.edu/sbia/software-agreement.html or COPYING file.

%  Author:
%  Erdem Varol
%  software@cbica.upenn.edu
%
%
% Reference: Varol, Erdem, Aristeidis Sotiras, Christos Davatzikos. 
% "MIDAS: regionally linear multivariate discriminative statistical mapping." NeuroImage (2018)

function map = midas(varargin)

if nargin==0
    printhelp()
    return
end

if( strcmp(varargin{1},'--help') || isempty(varargin))
    printhelp()
    return;
end

if( strcmp(varargin{1},'-h') || isempty(varargin) )
    printhelp()
    return
end

if( strcmp(varargin{1},'--version') || isempty(varargin) )
    fprintf('Version 1.0.\n')
    return
end

if( strcmp(varargin{1},'-v') || isempty(varargin) )
    fprintf('Version 1.0.\n')
    return
end

if( strcmp(varargin{1},'-u') || isempty(varargin) )
    fprintf(' EXAMPLE USE (in matlab) \n');
    fprintf(' midas(''-i'',''test.csv'',''-o'',''.'',''-r'',15,''-p'',200,''-c'',0.1) \n');
    fprintf(' EXAMPLE USE (in command line) \n');
    fprintf(' midas -i test.csv -o . -r 15 -p 200 -c 0.1 \n');
    return
end

if( strcmp(varargin{1},'--usage') || isempty(varargin) )
    fprintf(' EXAMPLE USE (in matlab) \n');
    fprintf(' midas(''-i'',''test.csv'',''-o'',''.'',''-r'',15,''-p'',200,''-c'',0.1) \n');
    fprintf(' EXAMPLE USE (in command line) \n');
    fprintf(' midas -i test.csv -o . -r 15 -p 200 -c 0.1 \n');
    return
end


% function returns MIDAS statistical maps and associated p-values computed
% by analytically estimating permutation testing
%
% INPUT
%
% REQUIRED
% [--input, -i] : .csv file containing full paths to input images. (REQUIRED)
%              We assume that the first column contains subject identifying
%              information; the second column contains the path to the
%              images, while the last column contains label information.
%              First line of the file should contain header information.

% [--outputDir, -o] : directory where the output from all folds will be saved (REQUIRED)
%
% OPTIONAL
%
%
% [--c, -c] : regularization parameter (positive scalar) (default C=0.1)
% [--radius, -r] : neighborhood radius in voxels (positive scalar) (default R=15)
% [--num, -p] : number of neighborhoods (positive scalar) (default P = 200 )
% [--usage, -u]  Prints basic usage message.          
% [--help, -h]  Prints help information.
% [--version, -v]  Prints information about software version.

%
% OUTPUT:
% map = structure that stores MIDAS statistics and p-values for every group/covariate
% provided with the following structure
% map.stat = cell that stores MIDAS statistic for groups/covariate in same
% order (output as .nii.gz file as well)
% map.p = cell that stores MIDAS statistic for groups/covariate in same
% order (output as .nii.gz file as well)
% map.N = neighborhood coverage amounts (output as .nii.gz file as well)
%
%
% NOTE: to compile this function do
% mcc -m  midas -A [NIFTI toolbox directory]
% EXAMPLE USE (in matlab)
% midas('-i','test.csv','-o','.','-r',15,'-p',200,'-c',0.1)
% EXAMPLE USE (in command line)
% midas -i test.csv -o . -r 15 -p 200 -c 0.1


if( sum(or(strcmpi(varargin,'--input'),strcmpi(varargin,'-i')))==1)
    niiCSV=varargin{find(or(strcmpi(varargin,'--input'),strcmp(varargin,'-i')))+1};
else
    error('run_midas_experiment_csv:argChk','Please specify input csv file!');
end


if( sum(or(strcmpi(varargin,'--outputDir'),strcmpi(varargin,'-o')))==1)
    outputDir=varargin{find(or(strcmp(varargin,'--outputDir'),strcmp(varargin,'-o')))+1};
else
    error('run_midas_experiment_csv:argChk','Please specify output directory!');
end

if( sum(or(strcmpi(varargin,'--c'),strcmpi(varargin,'-c')))==1)
    params.ls_svm_C=varargin{find(or(strcmpi(varargin,'--c'),strcmp(varargin,'-c')))+1};
else
    params.ls_svm_C=0.1;
end

if( sum(or(strcmpi(varargin,'--radius'),strcmpi(varargin,'-r')))==1)
    params.radius=varargin{find(or(strcmpi(varargin,'--radius'),strcmp(varargin,'-r')))+1};
else
    params.radius=15;
end

if( sum(or(strcmpi(varargin,'--num'),strcmpi(varargin,'-p')))==1)
    params.num_neighborhoods=varargin{find(or(strcmpi(varargin,'--num'),strcmp(varargin,'-p')))+1};
else
    params.num_neighborhoods=200;
end

if( sum(or(strcmpi(varargin,'--verbose'),strcmpi(varargin,'-vo')))==1)
    params.vo=varargin{find(or(strcmpi(varargin,'--verbose'),strcmp(varargin,'-vo')))+1};
else
    params.vo=0;
end

% create output directory
if (~exist(outputDir,'dir'))
    [status,~,~] = mkdir(outputDir);
    if (status == 0)
        error('run_midas_experiment_csv:argChk','Cannot create output directory!');
    end
end

params.ls_svm_C=input2num(params.ls_svm_C);
params.radius=input2num(params.radius);
params.num_neighborhoods=input2num(params.num_neighborhoods);
params.vo=input2num(params.vo);


% confirm validity of optional input arguments
validateFcn_C = @(x) (x>0);
validateFcn_R = @(x) (x>0);
validateFcn_P = @(x) isscalar(x) && (x>0) && (mod(x,1)==0);
validateFcn_vo = @(x) (x==0) || (x == 1);

if(~validateFcn_C(params.ls_svm_C))
    error('run_midas_experiment_csv:argChk','Regularization parameter should be positive!');
end

if(~validateFcn_R(params.radius))
    error('run_midas_experiment_csv:argChk','Radius parameter should be positive!');
end

if(~validateFcn_P(params.num_neighborhoods))
    error('run_midas_experiment_csv:argChk','Number of neighborhoods should be a positive integer!');
end

if(~validateFcn_vo(params.vo))
    error('run_midas_experiment_csv:argChk','VO parameter should be either 0 or 1!');
end

disp('Done');
disp('MIDAS runs with the following parameteres');
disp(['niiCSV: ' niiCSV]);
disp(['OutputDir: ' outputDir]);
disp(['C: ' num2str(params.ls_svm_C)]);
disp(['R: ' num2str(params.radius)]);
disp(['P: ' num2str(params.num_neighborhoods)]);
disp(['vo: ' num2str(params.vo)]);

% csv with features
fname=niiCSV;
if (~exist(fname,'file'))
    error('run_midas_experiment_csv:argChk','Input feature .csv file does not exist');
end


% input data
% assumption is that the first column contains IDs, the second paths to
% files and the last contains labels
disp('Loading feature images...');
input=readtable(fname);
ID=input{:,1};
fnames=input{:,2};
Y=input{:,3:end};
covariates=input.Properties.VariableNames(3:end);

% read the images
count = size(Y,1);
info = load_untouch_nii_gz(fnames{1});

img.dimx = info.hdr.dime.dim(2) ;
img.dimy = info.hdr.dime.dim(3) ;
img.dimz = info.hdr.dime.dim(4) ;

X = zeros(count,img.dimx*img.dimy*img.dimz);
template = load_untouch_nii_gz(fnames{1});
foreground=zeros(size(template));
for ii=1:count
    disp(['Loading image ' fnames{ii} ]);
    nii = load_untouch_nii_gz(fnames{ii});
    X(ii,:) = nii.img(:) ;
    foreground=foreground+abs(nii.img);
end
foreground=double(foreground>0);


% z-score imaging features
X=zscore(X);

for i=1:size(X,1)
    i
    data{i}=reshape(X(i,:),size(foreground));
end

clear X

toremove=find(isnan(sum(Y,2)));
if ~isempty(toremove)
    disp('Removed the following subjects due to incomplete information:')
    disp(ID(toremove));
    data(isnan(sum(Y,2)))=[];
    Y(isnan(sum(Y,2)),:)=[];
end
rng('shuffle')
map=midas_solver(data,Y,params,foreground);


disp('Saving results...')
if(params.vo==0)
    save([outputDir '/MIDAS_results.mat'],'map');
    for i=1:length(covariates)
        mat2nii(map.stat{1},fnames{1},[outputDir '/midas_stat_map_' covariates{i}]);
        mat2nii(map.p{1},fnames{1},[outputDir '/midas_p_map_' covariates{i}]);
    end
    mat2nii(map.N,fnames{1},[outputDir '/midas_coverage_map']);
else
    save([outputDir '/MIDAS_results.mat'],'map');
    for i=1:length(covariates)
        mat2nii(map.stat{1},fnames{1},[outputDir '/midas_stat_map_' covariates{i}]);
        mat2nii(map.p{1},fnames{1},[outputDir '/midas_p_map_' covariates{i}]);
    end
    mat2nii(map.N,fnames{1},[outputDir '/midas_coverage_map']);
end
disp('Done')
end

function mat2nii(map,template,prefix)
%function to write matlab map to a nifti file
%Input:
% 1) map: matlab map
% 2)template: a file to borrow nifti header from, one that matches the map
% dimensions
% 3) prefix: prefix of output
%Output:
% prefix.nii
% Example usage:
%mat2nii(map{1}.stat{1},'/cbica/software/external/fsl/4.1.5/data/standard/MNI152lin_T1_2mm_brain.nii.gz','output')

tmp=load_untouch_nii_gz(template);
tmp.hdr.dime.datatype=16;
tmp.img=map;

save_untouch_nii(tmp,prefix);
gzip([prefix '.nii'])
delete([prefix '.nii'])
disp(['Saved as ' prefix '.nii.gz'])
end

function printhelp()

fprintf(' function returns MIDAS statistical maps and associated p-values computed\n');
fprintf(' by analytically estimating permutation testing\n');
fprintf('\n');
fprintf(' INPUT\n');
fprintf('\n');
fprintf(' REQUIRED\n');
fprintf(' [--input, -i] : .csv file containing full paths to input images. (REQUIRED)\n');
fprintf('              We assume that the first column contains subject identifying\n');
fprintf('              information; the second column contains the path to the\n');
fprintf('              images, while the last column contains label information.\n');
fprintf('              First line of the file should contain header information.\n');
fprintf('\n');
fprintf(' [--outputDir, -o] : directory where the output from all folds will be saved (REQUIRED)\n');
fprintf('\n');
fprintf(' OPTIONAL\n');
fprintf('\n');
fprintf('\n');
fprintf(' [--c, -c] : regularization parameter (positive scalar) (default C=0.1)\n');
fprintf(' [--radius, -r] : neighborhood radius in voxels (positive scalar) (default R=15)\n');
fprintf(' [--num, -p] : number of neighborhoods (positive scalar) (default P = 200 )\n');
fprintf(' [--usage, -u]  Prints basic usage message.     \n');     
fprintf(' [--help, -h]  Prints help information.\n');
fprintf(' [--version, -v]  Prints information about software version.\n');
fprintf('\n');
fprintf('\n');
fprintf(' OUTPUT:\n');
fprintf(' map = structure that stores MIDAS statistics and p-values for every group/covariate\n');
fprintf(' provided with the following structure\n');
fprintf(' map.stat = cell that stores MIDAS statistic for groups/covariate in same\n');
fprintf(' order (output as .nii.gz file as well)\n');
fprintf(' map.p = cell that stores MIDAS statistic for groups/covariate in same\n');
fprintf(' order (output as .nii.gz file as well)\n');
fprintf(' map.N = neighborhood coverage amounts (output as .nii.gz file as well)\n');
fprintf('\n');
fprintf('\n');
fprintf(' NOTE: to compile this function do\n');
fprintf(' mcc -m  midas -A [NIFTI toolbox directory]\n');
fprintf(' EXAMPLE USE (in matlab)\n');
fprintf(' midas(''-i'',''test.csv'',''-o'',''.'',''-r'',15,''-p'',200,''-c'',0.1)\n');
fprintf(' EXAMPLE USE (in command line)\n');
fprintf(' midas -i test.csv -o . -r 15 -p 200 -c 0.1\n');
fprintf('======================================================================\n');
fprintf('Contact: software@cbica.upenn.edu\n');
fprintf('\n');
fprintf('Copyright (c) 2018 University of Pennsylvania. All rights reserved.\n');
fprintf('See COPYING file or http://www.med.upenn.edu/sbia/software/license.html\n');
fprintf('======================================================================\n');
end

function o=input2num(x)
if isnumeric(x)
    o=x;
else
    o = str2double(x);
end
end
