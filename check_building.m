%% check_building.m
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
%% carnotUIBK version 2.0
% Copyright (c) 2016-2019, University of Innsbruck, Unit for Energy 
% Efficient Building.
%
% Author    Date         Description
% DS,EL     2017-03-12   initial revision v1.0
% DS,EL	    2019-01-24   updates for GUI v2.0

%% update thermalzone
try
    building.thermalzone(building.variant_thermalzone) = building.thermalzone(building.variant_thermalzone).update_thermalzone(building.geometry(building.variant_geometry),building,building.variant_gains);
catch ME
    ME
    err = ('An error occured during the update of the thermalzones for the simulation.');
    error(err)
    return
end

%% check building
% check zones
b2t = building.get_building();

for jj = 1:NUMBEROFZONES
    % ZONES
    if length(b2t.thermalzone.zone(jj).name) > 4 && strcmp(b2t.thermalzone.zone(jj).name(1:5), 'EMPTY')
        if ~isempty(building.get_building().thermalzone.zone(jj).matrix_wd) || ~isempty(building.get_building().thermalzone.zone(jj).matrix_wi) || ~isempty(building.get_building().thermalzone.zone(jj).matrix_gains)
            err = ['The zone ' num2str(jj) ' should be empty, but it seems not to be the case. Check the name of the zone. EMTPY is not allowed as a key-word.'];
            error(err)
            return
        end
    else
        for jk = 1:size(b2t.thermalzone.zone(jj).matrix_wd,1)
            if isempty(b2t.construction.get_structure(b2t.thermalzone.zone(jj).matrix_wd{jk,1}))
                err = ['The construction ''' b2t.thermalzone.zone(jj).matrix_wd{jk,1} ''' is not defined! Please define it in the variant ' num2str(building.variant_construction) '.'];
                error(err)
                return
            end
            if ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'AMBIENT') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'INTERNAL') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'GROUND') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'NEIGHBOUR') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'GROUND_1') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'NEIGHBOUR_1') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'GROUND_2') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'NEIGHBOUR_2') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'GROUND_3') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'NEIGHBOUR_3') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'GROUND_4') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'NEIGHBOUR_4') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'GROUND_5') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'NEIGHBOUR_5') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'GROUND_6') && ~strcmp(b2t.thermalzone.zone(jj).matrix_wd{jk,3},'NEIGHBOUR_6')
                err = ['The boundary condition ' b2t.thermalzone.zone(jj).matrix_wd{jk,3} ' in zone ' num2str(jj) ' is not well defined! Is the room existing or is there a typing error?'];
                error(err)
                return
            end
            if b2t.thermalzone.zone(jj).matrix_wd{jk,2} == 0.0
                err = ['The area of a wall in the zone ' num2str(jj) ' is zero. This is not allowed!'];
                error(err)
                return
            end
        end
        for jk = 1:size(b2t.thermalzone.zone(jj).matrix_wi,1)
            if isempty(b2t.construction.get_structure(b2t.thermalzone.zone(jj).matrix_wi{jk,1}))
                err = ['The construction ''' b2t.thermalzone.zone(jj).matrix_wi{jk,1} ''' is not defined! Please define it in the variant ' num2str(building.variant_construction) '.'];
                error(err)
                return
            end
            if ~strcmp(b2t.thermalzone.zone(jj).matrix_wi{jk,4},'AMBIENT')
                err = ['A window in the zone ' num2str(jj) ' is not well defined. It has not the boundary condition AMBIENT, but ' b2t.thermalzone.zone(jj).matrix_wi{jk,4} '.'];
                error(err)
                return
            end
            if b2t.thermalzone.zone(jj).matrix_wi{jk,2} == 0.0
                err = ['The area of a window in the zone ' num2str(jj) ' is zero. This is not allowed!'];
                error(err)
                return
            end
        end
    end

    % INTERSECTIONS
    if jj < 10
        for jjj = 1:NUMBEROFZONES
            if ~isempty(b2t.thermalzone.intersection{jj,jjj})
                if length(b2t.thermalzone.intersection{jj,jjj}.name) > 4 && strcmp(b2t.thermalzone.intersection{jj,jjj}.name(1:5), 'EMPTY')
                    if ~isempty(building.get_building().thermalzone.intersection{jj,jjj}.matrix_wd)
                        err = ['The intersection (' num2str(jj) ',' num2str(jjj) ') should be empty, but it seems not to be the case. Check the name of the zone. EMTPY is not allowed as a key-word.'];
                        error(err)
                        return
                    end
                else
                    for jk = 1:size(b2t.thermalzone.intersection{jj,jjj}.matrix_wd,1)
                        if isempty(b2t.construction.get_structure(b2t.thermalzone.intersection{jj,jjj}.matrix_wd{jk,1}))
                            err = ['The construction ''' b2t.thermalzone.intersection{jj,jjj}.matrix_wd{jk,1} ''' is not defined! Please define it in the variant ' num2str(building.variant_construction) '.'];
                            error(err)
                            return
                        end
                        if b2t.thermalzone.intersection{jj,jjj}.matrix_wd{jk,2} == 0.0
                            err = ['The area of a wall in the intersection (' num2str(jj) ',' num2str(jjj) ') is zero. This is not allowed!'];
                            error(err)
                            return
                        end
                    end
                end
            end
        end
    end
end

% building solar ratio
for jj = 1:NUMBEROFZONES
    if length(b2t.thermalzone.zone(jj).name) > 4 && strcmp(b2t.thermalzone.zone(jj).name(1:5), 'EMPTY')
        % zone empty!
    else
        b2t.thermalzone.check_solar_ratio(jj)
    end
end

%% plot matrizes to check
if PLOT
    % zones
    for jj = 1:NUMBEROFZONES
        b2t.thermalzone.check_matrix(jj)
    end

    % intersections
    for jj = 1:NUMBEROFZONES
        for jk = (jj+1):NUMBEROFZONES
            b2t.thermalzone.check_matrix(jj,jk)
        end
    end
end

%%
clear b2t