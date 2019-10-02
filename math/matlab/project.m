function [projectroot, figpath, varpath, analysispath] = project(project_name, varargin)
%PROJECT    Configure the current environment for the given project.
%

projectroot = [getenv('HOME') '/research/prj/' project_name];
varpath = [getenv('HOME') '/research/var/' project_name];
docpath = '';
analysispath = [projectroot '/analysis'];

if(nargin > 1);
    iend = size(varargin,2);
    i = 1;
    while(i<=iend)
        if(strcmp(varargin{i}, 'docpath'))
            docpath = ['/doc/' varargin{i+1}];
            i=i+1;
        elseif(strcmp(varargin{i}, 'projectroot'))
            projectroot = varargin{i+1};
            i=i+1;           
        elseif(strcmp(varargin{i}, 'varpath'))
            varpath = varargin{i+1};
            i=i+1;           
        end
        i=i+1;
    end
end

addpath([projectroot '/math'], '-begin');
figpath = [projectroot docpath '/figures/'];
cd(varpath);