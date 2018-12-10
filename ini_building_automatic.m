%% ini_building_automatic.m
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
function building = ini_building_automatic(building, modify, import_mode, name_XML, name_EXCEL, name_PHPP, NUMBEROFZONES)
    %% gbXML file (DO NOT MODIFY)
    building.XML.name = name_XML;
    
    %% EXCEL file (DO NOT MODIFY)
    building.EXCEL.name = name_EXCEL;
    
    %% PHPP file (DO NOT MODIFY)
    building.PHPP.name = name_PHPP;
    
    %% set variant numbers (DO NOT MODIFY)
    disp(' ')
    switch modify.geometry
        case -1
            disp('No new or updated variant for geometry.')
        case 0
            disp('New automatic variant for geometry.')
            if size(building.geometry(1).room,2) == 0
                building.variant_geometry = 1;
            else
                building.variant_geometry = length(building.geometry) + 1;
            end
            disp(['   -> new variant geometry: ' num2str(building.variant_geometry)])
        otherwise
            disp(['Overwrite variant ' num2str(modify.geometry) ' for geometry.'])
            building.variant_geometry = modify.geometry;
    end
    
    switch modify.construction
        case -1
            disp('No new or updated variant for construction.')
        case 0
            disp('New automatic variant for construction.')
            if size(building.geometry(1).room,2) == 0
                building.variant_construction = 1;
            else
                building.variant_construction = length(building.construction) + 1;
            end
            disp(['   -> new variant construction: ' num2str(building.variant_construction)])
        otherwise
            disp(['Overwrite variant ' num2str(modify.construction) ' for construction.'])
            building.variant_construction = modify.construction;
    end
    
    switch modify.boundary
        case -1
            disp('No new or updated variant for boundary.')
        case 0
            disp('New automatic variant for boundary.')
            if size(building.geometry(1).room,2) == 0
                building.variant_boundary = 1;
            else
                building.variant_boundary = length(building.boundary) + 1;
            end
            disp(['   -> new variant boundary: ' num2str(building.variant_geometry)])
        otherwise
            disp(['Overwrite variant ' num2str(modify.boundary) ' for boundary.'])
            building.variant_boundary = modify.boundary;
    end
    
    switch modify.gains
        case -1
            disp('No new or updated variant for gains.')
        case 0
            disp('New automatic variant for gains.')
            if size(building.geometry(1).room,2) == 0
                building.variant_gains = 1;
            else
                building.variant_gains = length(building.gains) + 1;
            end
            disp(['   -> new variant gains: ' num2str(building.variant_gains)])
        otherwise
            disp(['Overwrite variant ' num2str(modify.gains) ' for gains.'])
            building.variant_gains = modify.gains;
    end
    
    switch modify.hvac
        case -1
            disp('No new or updated variant for hvac.')
        case 0
            disp('New automatic variant for hvac.')
            if size(building.geometry(1).room,2) == 0
                building.variant_hvac = 1;
            else
                building.variant_hvac = length(building.hvac) + 1;
            end
            disp(['   -> new variant hvac: ' num2str(building.variant_hvac)])
        otherwise
            disp(['Overwrite variant ' num2str(modify.hvac) ' for hvac.'])
            building.variant_hvac = modify.hvac;
    end
    
    switch modify.thermalzone
        case -1
            disp('No new or updated variant for thermalzone.')
        case 0
            disp('New automatic variant for thermalzone.')
            if size(building.geometry(1).room,2) == 0
                building.variant_thermalzone = 1;
            else
                building.variant_thermalzone = length(building.thermalzone) + 1;
            end
            disp(['   -> new variant thermalzone: ' num2str(building.variant_thermalzone)])
        otherwise
            disp(['Overwrite variant ' num2str(modify.thermalzone) ' for thermalzone.'])
            building.variant_thermalzone = modify.thermalzone;
    end
    
    %% for PHPP: load parameters
    if strcmp(import_mode,'PHPP')
        [PHPP_parameter_num, ~, PHPP_parameter] = xlsread(building.PHPP.name,'CARNOT','C3:C12');
        
        % language of the PHPP (1: German, 0: English)
        building.PHPP.language = PHPP_parameter_num(1);
        
        % version of the PHPP (implemented for '9.1' version)
        building.PHPP.version = PHPP_parameter{2};
        
        % model of the construction of the wall: 'UA' or 'RC' (for model 'RC' in excel add 2 columns
        % with rho and c between M and N in the sheet "U-Werte" or "U-Values")
        building.PHPP.choice_modelconswall = PHPP_parameter{3};
        
        % choose the max heating power (as the PHPP: 'limited', or 'unlimited')
        building.PHPP.choice_heatload = PHPP_parameter{4};
        
        % weather
        building.PHPP.weather.name = PHPP_parameter{5};
        building.PHPP.weather.file = PHPP_parameter{6};
        building.PHPP.weather.latitude = PHPP_parameter_num(7);
        building.PHPP.weather.longitude = PHPP_parameter_num(8);
        building.PHPP.weather.ref_median = PHPP_parameter_num(9);
        
        % orientation important
        building.PHPP.orientation_important = PHPP_parameter_num(10);
    end
    
    %% define building (DO NOT MODIFY)
    % add the geometry (DO NOT MODIFY)
    if modify.geometry >= 0
        if strcmp(import_mode,'gbXML')
            excel_create = 'structure';
            building.geometry(building.variant_geometry) = GEOMETRY();
            building = building.add_variant_geometry(building.variant_geometry, building.geometry(building.variant_geometry).geometry_from_XML_andor_excel(building.XML.name, building, building.EXCEL.name, excel_create));
        elseif strcmp(import_mode,'excel')
            building.geometry(building.variant_geometry) = GEOMETRY();
            building = building.add_variant_geometry(building.variant_geometry, building.geometry(building.variant_geometry).geometry_from_excel(building.EXCEL.name, building));
        elseif strcmp(import_mode,'PHPP')
            building.geometry(building.variant_geometry) = GEOMETRY();
            building = building.add_variant_geometry(building.variant_geometry, building.geometry(building.variant_geometry).geometry_from_PHPP(building.PHPP.name, building.PHPP.language, building.PHPP.version, building.PHPP.choice_modelconswall, building));
            % 13-11-2018 Eleonora Leonardi
            building.geometry.geometry_to_excel(['EXCEL_' building.PHPP.name], building, building.variant_geometry);
        end
    end
    
    % add the construction (MODIFY the function "create_construction")
    if modify.construction >= 0
        if strcmp(import_mode,'gbXML')
            cons = CONSTRUCTION();
            cons = cons.construction_from_EXCEL(building.EXCEL.name);
            cons = cons.construction_windows_from_EXCEL(building.EXCEL.name);
            building = building.add_variant_construction(building.variant_construction, cons);
        elseif strcmp(import_mode,'excel')
            cons = CONSTRUCTION();
            cons = cons.construction_from_EXCEL(building.EXCEL.name);
            cons = cons.construction_windows_from_EXCEL(building.EXCEL.name);
            building = building.add_variant_construction(building.variant_construction, cons);
        elseif strcmp(import_mode,'PHPP')
            cons = CONSTRUCTION();
            cons = cons.construction_from_PHPP(building.PHPP.name, building.PHPP.language, building.PHPP.version, building.PHPP.choice_modelconswall);
            cons = cons.constructionwindow_from_PHPP(building.PHPP.name, building.PHPP.language, building.PHPP.version);
            cons = cons.constructiontbUA_from_PHPP(building.PHPP.name, building.PHPP.language, building.PHPP.version);
            building = building.add_variant_construction(building.variant_construction, cons);
            % 13-11-2018 Eleonora Leonardi
            building.construction.construction_to_excel(['EXCEL_' building.PHPP.name], building, building.variant_construction);
        end
    end
    
    % add the gains
    if modify.gains >= 0
        if strcmp(import_mode,'gbXML')
            intgain = GAINS();
            building = building.add_variant_gains(building.variant_gains, intgain.gains_from_excel(building.EXCEL.name));
        elseif strcmp(import_mode,'excel')
            intgain = GAINS();
            building = building.add_variant_gains(building.variant_gains, intgain.gains_from_excel(building.EXCEL.name));
        elseif strcmp(import_mode,'PHPP')
            intgain = GAINS();
            intgain = intgain.gain_from_PHPP(building.PHPP.name, building.PHPP.language, building.PHPP.version);
            intgain = intgain.assign_gain_to_zone('Person_W/m²_from_PHPP', 1);
            building = building.add_variant_gains(building.variant_gains, intgain);
            % 13-11-2018 Eleonora Leonardi
            building.gains.gains_to_excel(['EXCEL_' building.PHPP.name], building, building.variant_gains);
        end
    end
    
    % add the thermal zone
    if modify.thermalzone >= 0
        if strcmp(import_mode,'gbXML')
            thzo = THERMALZONE(NUMBEROFZONES);
            building = building.add_variant_thermalzone(building.variant_thermalzone, thzo.zone_from_XML(building.XML.name));
        elseif strcmp(import_mode,'excel')
            thzo = THERMALZONE(NUMBEROFZONES);
            building = building.add_variant_thermalzone(building.variant_thermalzone, thzo.zone_fromexcel(building.EXCEL.name));
        elseif strcmp(import_mode,'PHPP')
            thzo = THERMALZONE(NUMBEROFZONES);
            if strcmp(building.PHPP.choice_modelconswall, 'UA')
                model_room = 2;
            elseif strcmp(building.PHPP.choice_modelconswall, 'RC')
                model_room = 3;
            end
            thzo = thzo.add_zone(1, 'building', {'whole_building'}, model_room);
            thzo = thzo.add_cpspectozone_setorientimp(building.PHPP.name, building.PHPP.language, building.PHPP.version, 1, building.PHPP.orientation_important);
            building = building.add_variant_thermalzone(building.variant_thermalzone, thzo);
            % 13-11-2018 Eleonora Leonardi
            building.thermalzone.zone_to_excel(['EXCEL_' building.PHPP.name], building, building.variant_thermalzone);
        end
    end
    
    % add the boundary conditions
    if modify.boundary >= 0
        if strcmp(import_mode,'gbXML')
            boun = BOUNDARY();
            boun = boun.add_ground_neighbour_weather_from_excel(building, building.EXCEL.name);
            boun = boun.complete_ground_and_neighbour(building);
            building = building.add_variant_boundary(building.variant_boundary, boun);
        elseif strcmp(import_mode,'excel')
            boun = BOUNDARY();
            boun = boun.add_ground_neighbour_weather_from_excel(building, building.EXCEL.name);
            boun = boun.complete_ground_and_neighbour(building);
            building = building.add_variant_boundary(building.variant_boundary, boun);
        elseif strcmp(import_mode,'PHPP')
            boun = BOUNDARY();
            boun = boun.neigbour_4_from_PHPP(building, building.PHPP.name, building.PHPP.language, building.PHPP.version);
            boun = boun.ground_temp_from_PHPP(building, building.PHPP.name, building.PHPP.language, building.PHPP.version);
            boun = boun.add_weather_from_file(building, building.PHPP.weather.file, building.PHPP.weather.name, building.PHPP.weather.latitude, building.PHPP.weather.longitude, building.PHPP.weather.ref_median);
            boun = boun.complete_ground_and_neighbour(building);
            building = building.add_variant_boundary(building.variant_boundary, boun);
            % 13-11-2018 Eleonora Leonardi
            building.boundary.ground_neighbour_weather_to_excel(['EXCEL_' building.PHPP.name], building, building.variant_boundary);
        end
    end
    
    % add an empty hvac (except for PHPP)
    if modify.hvac >= 0
        if strcmp(import_mode,'gbXML')
            hvac = HVAC();
            building = building.add_variant_hvac(building.variant_hvac, hvac);
        elseif strcmp(import_mode,'excel')
            hvac = HVAC();
            building = building.add_variant_hvac(building.variant_hvac, hvac);
        elseif strcmp(import_mode,'PHPP')
            hvac = HVAC();
            hvac = hvac.add_system_from_PHPP('HVAC from PHPP', 1, building.PHPP.name, building.PHPP.language, building.PHPP.version, building.PHPP.choice_heatload);
            building = building.add_variant_hvac(building.variant_hvac, hvac);
        end
    end
end
