%% load_building_automatic.m
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
%% carnotUIBK version 1.3
% Copyright (c) 2016-2018, University of Innsbruck, Unit for Energy 
% Efficient Building.
%
% Author    Date         Description
% DS        2018-05-25   v1.3: creation of a more simplified way to load and modify the building
%%
function [PLOT NUMBEROFZONES NUMBEROFINTERSECTIONS NUMBEROFWALLSINZONES NUMBEROFWINDOWSINZONES NUMBEROFWALLSININTERSECTIONS NUMBEROFWINDOWSINZONESANDINTERSECTIONS NUMBEROFGAINSINZONES NUMBEROFCONTROLS building] = load_building_automatic(name_build_obj)
    %% emtpy building object and global variables (DO NOT MODIFY)
    building = BUILDING();

    PLOT = 0;
    NUMBEROFZONES = 10;
    NUMBEROFINTERSECTIONS = (NUMBEROFZONES*NUMBEROFZONES-NUMBEROFZONES)/2;
    NUMBEROFWALLSINZONES = 10;
    NUMBEROFWINDOWSINZONES = 5;
    NUMBEROFWALLSININTERSECTIONS = 3;
    NUMBEROFWINDOWSINZONESANDINTERSECTIONS = NUMBEROFWALLSINZONES+(NUMBEROFZONES-1)*NUMBEROFWALLSININTERSECTIONS;
    NUMBEROFGAINSINZONES = 5;
    NUMBEROFCONTROLS = 20;

    %% define building name (DO NOT MODIFY)
    building.name = name_build_obj;

    %% load or create building (DO NOT MODIFY)
    if exist([building.name '.mat'] , 'file') == 2
        load(building.name);
    else
        building = building.save();
    end
end