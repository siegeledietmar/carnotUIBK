%% init_carnotUIBK.m
% ***********************************************************************
% This file is part of the uibkCARNOT Blockset.
% 
% Copyright (c) 2016-2019, University of Innsbruck, Unit for Energy 
% Efficient Building.
%   Dietmar Siegele     dietmar.siegele@uibk.ac.at
%   Eleonora Leonardi   eleonora.leonardi@uibk.ac.at
% Additional Copyright for this file see list auf authors.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice, 
%    this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright 
%    notice, this list of conditions and the following disclaimer in the 
%    documentation and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its 
%    contributors may be used to endorse or promote products derived from 
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
% THE POSSIBILITY OF SUCH DAMAGE.
% **********************************************************************
% 
%% carnotUIBK
% Copyright (c) 2016-2019, University of Innsbruck, Unit for Energy 
% Efficient Building.

% DS,EL	    2019-01-24   updates for GUI v2.0

%% this ist a copy of the CARNOT installation
fprintf('################################################\n*Initializing carnotUIBK Toolbox\n\n')

cpath = pwd;
ctrl = 'setpaths';

% find root path of carnot
if exist('carnotUIBK','file') == 4
    rootpath = fileparts(which('carnotUIBK'));
else % paths might not be added yet, so carnotUIBK.slx won't be found
    warning('MATLAB:Path', ...
        'carnotUIBK.slx could not be found. Possibly its path has not been added yet. Trying to continue with alternative method.')
    ptemp = pwd;
    cd(fileparts(mfilename('fullpath')));
    cd('..\..') % supposed this script is located in <root>\public\common\scripts--> changed to <root>\public\src_m
    if exist('carnotUIBK','file') ~= 4
        error('carnotUIBK.slx could not be found. Check consistency of subpaths')
    end
    rootpath = pwd;
    cd(ptemp)
end

% define standardized carnotUIBK paths
carnotUIBKpaths = {...
    'root'          fullfile(rootpath);...
    'TEMPLATE'      fullfile(rootpath,'TEMPLATE');...
    'result_scripts'      fullfile(rootpath,'result_scripts');...
    'HVAC_SYSTEM'      fullfile(rootpath,'HVAC_SYSTEM');...
    };

% set paths mode: adds carnotUIBK paths to matlab path
if strcmp(ctrl,'setpaths')
    % specify which paths should be added
    paths2add = { ...               % first to add is last in path
        %'help','inthelp', ...    
        'root', ...
        'result_scripts', ...
        'HVAC_SYSTEM', ...
        };
    disp('Adding carnotUIBK paths...')
    for i = 1:length(paths2add)
        for j = 1:size(carnotUIBKpaths,1)
            if strcmp(carnotUIBKpaths(j,1),paths2add(i))
                addpath(cell2mat(carnotUIBKpaths(j,2)))
                disp(['  ... ' cell2mat(carnotUIBKpaths(j,2))])
                break
            end
        end
    end
else % return demanded path as string
    for i = 1:size(carnotUIBKpaths,1)
        if strcmp(carnotUIBKpaths(i,1),ctrl)
            p = cell2mat(carnotUIBKpaths(i,2));
            break
        end
    end
end

rehash;

cd(cpath)

fprintf('\n*done\n\nRead the README for further instructions.\n################################################\n')

%% save PATH
disp('save path')
savepath;

%% update Toolbox cache
disp('update toolbox cache')
rehash toolboxcache;

cd(cpath)

%%
clear
