%% GENERATECODE.m
% ***********************************************************************
% This file is part of the uibkCARNOT Blockset.
% 
% Copyright (c) 2016-2018, University of Innsbruck, Unit for Energy 
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
%% carnotUIBK version 1.3.2
% Copyright (c) 2016-2018, University of Innsbruck, Unit for Energy 
% Efficient Building.
%
% Author    Date         Description
% DS,EL     2017-03-12   initial revision v1.0
% DS        2018-05-28   v1.3: remove 'wall_node_optim'
% DS        2018-07-02   v1.3.1: bugs in PHPP
% DS        2018-07-11   v1.3.2: bugs in PHPP

%% file used for compiling files and help libary
securefiles = 0;

% files needs to be compiled
files2compile = {'BUILDING','BOUNDARY','CONSTRUCTION','DOOR','GAIN','GAIN_TO_ZONE','GAINS','GEOMETRY','GROUND','HVAC','INTERSECTION','NEIGHBOUR','RES','ROOM','STRUCTURE','SYSTEM_HVAC','THERMALZONE','WALL','WEATHER','WINDOW','ZONE','create_building_automatic','ini_building_automatic','load_building_automatic'};

% check directories and create them
if ~exist('carnotUIBK', 'dir')
  mkdir('carnotUIBK');
end
if ~exist('carnotUIBK/html', 'dir')
  mkdir('carnotUIBK/html');
end

% compile m-files and copy
for jj = 1:length(files2compile)
    if securefiles
        pcode(files2compile{jj});
        movefile([files2compile{jj} '.p'],['carnotUIBK\' files2compile{jj} '.p']);
    else
        copyfile([files2compile{jj} '.m'],['carnotUIBK\' files2compile{jj} '.m']);
    end
end
pcode('run4simulink');
movefile(['run4simulink.p'],['carnotUIBK\run4simulink.p']);

% deploy html-help and copy
for jj = 1:length(files2compile)
    htmlstr = help2html(files2compile{jj},'','-doc');
    fid = fopen([files2compile{jj} '.html'],'w');
    fprintf(fid,'%s',htmlstr);
    fclose(fid);
    movefile([files2compile{jj} '.html'],['carnotUIBK\html\' files2compile{jj} '.html']);
end

% copy carnotUIBK simulink libary
copyfile('carnotUIBK.slx','carnotUIBK\carnotUIBK.slx')

% copy specific compiled files of libary
copyfile('pdepesim_hygro_c_v4.mexw64','carnotUIBK\pdepesim_hygro_c_v4.mexw64')
copyfile('cprintf.m','carnotUIBK\cprintf.m')
copyfile('xml2struct.m','carnotUIBK\xml2struct.m')
copyfile('xlsread1.m','carnotUIBK\xlsread1.m')
copyfile('xlswrite1.m','carnotUIBK\xlswrite1.m')

% copy sketchup
% copyfile('SKETCHUP','carnotUIBK\SKETCHUP')

% copy templates
copyfile('TEMPLATE','carnotUIBK\TEMPLATE')
% copyfile('template_XLS.xlsx','carnotUIBK\template_XLS.xlsx')
% copyfile('template_XML.xlsx','carnotUIBK\template_XML.xlsx')

% copy manual
copyfile('User Manual_v1.pdf','carnotUIBK\User Manual_v1.pdf')
