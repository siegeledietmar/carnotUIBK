%% BOUNDARY.m
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
% DS        2017-03-16   v1.1: fixed several bugs for empty cells reading
%                        PHPP file caused troubles with NaNs
% DS,EL	    2019-01-24   updates for GUI v2.0

%%
classdef BOUNDARY
    % BOUNDARY
    
    properties
        name = 'new boundary';
        weather = [];
        ground = GROUND();
        neighbour = NEIGHBOUR();
    end
    
    methods
        function obj = BOUNDARY()
            
        end
        
        function obj = add_neighbour_set_temperature(obj, building, number, name, model, temperature_time, temperature_value)
            % 1 ... building object
            % 2 ... name
            % 3 ... model -1 no model
                 % 0 fixed value of temperature for the whole year
                 % 1 one fixed value of temperature for the winter (first) and one fixed value of temperature for the summer (second)
                 % 2 sequence of values with time (give it for one year)
            % 4 ... temperature_time: 1 value of temperature (model 0) or 2 value of temperature, the first for the winter and the second for the summer (model 1) or time serie connected to value serie (model 2)
            % 5 ... optional: temperature_value: none (model 0 and 1) or value serie connected to time serie (model 2)
            if (nargin == 6)
                obj.neighbour(number) = NEIGHBOUR(building, name, model, temperature_time);

            elseif (nargin == 7)
                obj.neighbour(number) = NEIGHBOUR(building, name, model, temperature_time, temperature_value);
            else
                error('Not valid NEIGHBOUR, too much or too less parameters')
            end
        end
        
        function obj = add_ground_set_temperature(obj, building, number, name, model, temperature_time, temperature_value)
            % 1 ... building object
            % 2 ... name
            % 3 ... model -1 no model
                 % 0 fixed value of temperature for the whole year
                 % 1 one fixed value of temperature for the winter (first) and one fixed value of temperature for the summer (second)
                 % 2 sequence of values with time (give it for one year)
            % 4 ... temperature_time: 1 value of temperature (model 0) or 2 value of temperature, the first for the winter and the second for the summer (model 1) or time serie connected to value serie (model 2)
            % 5 ... temperature_value: none (model 0 and 1) or value serie connected to time serie (model 2)
            if (nargin == 6)
                obj.ground(number) = GROUND(building, name, model, temperature_time);
            elseif (nargin == 7)
                obj.ground(number) = GROUND(building, name, model, temperature_time, temperature_value);
            else
                error('Not valid GROUND, too much or too less parameters')
            end
        end
        
        function obj = ground_temp_from_PHPP(obj, building, name_xls_PHPP, language, version)
            % 1 ... building object
            % 2 ... name of the file xls of the PHPP
            % 3 ... language: 1:german, 0:english
            % 4 ... version of PHPP used
            
            filename = name_xls_PHPP;
            
            switch version
                case '9.1'
%                     if language
%                         sheet = 'Erdreich';
%                     else
%                         sheet = 'Ground';
%                     end
%                     range = 'E110:P111';
%                     [~, ~, data] = xlsread(filename, sheet, range);
%                     index_raw_win = 1;
%                     index_column_win = 1:12;
%                     index_raw_sum = 2;
%                     index_column_sum = 1:12;
%             end
%             for ii = 1:12
%             temp_ground_win(ii) = data{index_raw_win,index_column_win(ii)};
%             temp_ground_sum(ii) = data{index_raw_sum,index_column_sum(ii)};
%             end
%             temperature_time = [0 31 59 90 120 151 181 212 243 273 304 334]*24*3600;
%             temperature_value = [temp_ground_win(1:4) temp_ground_sum(5:9) temp_ground_win(10:12)];
%             obj = obj.add_ground_set_temperature(building, 1, 'ground_from_PHPP', 2, temperature_time, temperature_value);
                    if language
                        sheet = 'Klima';
                    else
                        sheet = 'Climate';
                    end
                    range = 'E32:P32';
                    [~, ~, data] = xlsread(filename, sheet, range);
            end
            for ii = 1:12
            temp_ground(ii) = data{1,ii};
            end
            temperature_time = [0 31 59 90 120 151 181 212 243 273 304 334]*24*3600;
            temperature_value = temp_ground;
            obj = obj.add_ground_set_temperature(building, 1, 'ground_from_PHPP', 2, temperature_time, temperature_value);
        end
        
        function obj = neigbour_4_from_PHPP(obj, building, name_xls_PHPP, language, version)
            % 1 ... building object
            % 2 ... name of the file xls of the PHPP
            % 3 ... language: 1:german, 0:english
            % 4 ... version of PHPP used
            
            filename = name_xls_PHPP;
            
             switch version
                case '9.1'
                    if language
                        sheet = 'Klima';
                    else
                        sheet = 'Climate';
                    end
                    range = 'E24:P24';
                    [~, ~, data] = xlsread(filename, sheet, range);
                    range = 'E32:P32';
                    [~, ~, data_B] = xlsread(filename, sheet, range);
                    
                    if language
                        sheet = 'Flächen';
                    else
                        sheet = 'Areas';
                    end
                    range = 'AB21';
                    [~, ~, data1] = xlsread(filename, sheet, range);
                    range = 'K19:K20';
                    [~, ~, data1_B] = xlsread(filename, sheet, range);
                    
                    if language
                        sheet = 'Nachweis';
                    else
                        sheet = 'Verification';
                    end
                    range = 'K27:N29';
                    [~, ~, data2] = xlsread(filename, sheet, range);
                    index_raw_setpointwin = 1;
                    index_column_setpointwin = 1;
                     index_raw_setpointsum = 1;
                     index_column_setpointsum = 4;
             end
            
             
            teta_i(1:12) = data2{index_raw_setpointwin,index_column_setpointwin};
            for ii = 1:12
                teta_e(ii) = data{1,ii};
            end
            red_fac = data1{1,1};
            delta_t(1:12) = (teta_i-teta_e)*red_fac;
            temperature_time = [0 31 59 90 120 151 181 212 243 273 304 334]*24*3600;
            temperature_value(1:12) = teta_i(1:12) - delta_t(1:12);
            if sum(isnan(temperature_value))
                warning('Neighbour temperature 4 set to 20°C')
                obj = obj.add_neighbour_set_temperature(building, 4, 'neighbour_with_red_factor', 2, temperature_time, ones(1,12)*20);
            else
                obj = obj.add_neighbour_set_temperature(building, 4, 'neighbour_with_red_factor', 2, temperature_time, temperature_value);
            end
            
            %definition ground temperature
            for ii = 1:12
                temp_ground(ii) = data_B{1,ii};
            end
            
            %definition internal temperature
            for jj=1:5
                temp_int(jj)= data2{index_raw_setpointwin,index_column_setpointwin};
            end
            
            for jj=6:8
                temp_int(jj)= data2{index_raw_setpointsum,index_column_setpointsum};
            end
            
            for jj=9:12
                temp_int(jj)= data2{index_raw_setpointwin,index_column_setpointwin};
            end
            
            % NEIGHBOUR_1 is I
            obj = obj.add_neighbour_set_temperature(building, 1, 'internal temperature', 2, temperature_time, temp_int);
            
            switch data1_B{1,1}
                case 'B'
                    obj = obj.add_neighbour_set_temperature(building, 2, 'ground temperature', 2, temperature_time, temp_ground);
                    
                case 'A'
                    obj = obj.add_neighbour_set_temperature(building, 2, 'external temperature', 2, temperature_time, teta_e);
                    
                case 'P'
                    obj = obj.add_neighbour_set_temperature(building, 2, 'ground temperature', 2, temperature_time, temp_ground);
                    
                case 'X'
                    obj = obj.add_neighbour_set_temperature(building, 2, 'neighbour_with_red_factor', 2, temperature_time, temperature_value);
            end
            
            switch data1_B{2,1}
                case 'B'
                    obj = obj.add_neighbour_set_temperature(building, 3, 'ground temperature', 2, temperature_time, temp_ground);
                    
                case 'A'
                    obj = obj.add_neighbour_set_temperature(building, 3, 'external temperature', 2, temperature_time, teta_e);
                    
                case 'P'
                    obj = obj.add_neighbour_set_temperature(building, 3, 'ground temperature', 2, temperature_time, temp_ground);
                    
                case 'X'
                    obj = obj.add_neighbour_set_temperature(building, 3, 'neighbour_with_red_factor', 2, temperature_time, temperature_value);
            end
            
        end
        
        function obj = add_ground_neighbour_weather_from_excel(obj, building, name_xls)
            % to take the boundary conditions from the excel (only template
            % excel supported)
            % 1 ... building object
            % 2 ... name of the excel file
            vollpfad = [pwd '\' name_xls];
            
            if exist(vollpfad, 'file')
                warning('Existing Excel file used.')
            else
                rootpath = fileparts(which('carnotUIBK'));
                copyfile([rootpath '\building_TEMPLATE.xlsx'], vollpfad);
                warning('New Excel file generated.')
            end
            
            [~, ~, raw_boundary] = xlsread(name_xls, 'Boundary');
            
            raw_ground = raw_boundary(4:9,:);
            raw_neighbour = raw_boundary(14:19,:);
            raw_weather = raw_boundary(23,:);
            
            % to add the ground temperature
            for ii = 1:size(raw_ground,1)
                if isnan(raw_ground{ii,2})
                    %break
                else
                    number = raw_ground{ii,1};
                    name = raw_ground{ii,2};
                    model = raw_ground{ii,3};
                    temper_value = [];
                    if model == 0
                        temper_value = raw_ground{ii,4};
                        obj = obj.add_ground_set_temperature(building, number, name, model, temper_value);
                    elseif model == 1
                        temper_value = [raw_ground{ii,4} raw_ground{ii,5}];
                        obj = obj.add_ground_set_temperature(building, number, name, model, temper_value);
                    elseif model == 2
                        temper_time = [0 31 59 90 120 151 181 212 243 273 304 334]*24*3600;
                        for ll = 1:12
                            temper_value(1,ll) = raw_ground{ii,ll+3};
                        end
                        obj = obj.add_ground_set_temperature(building, number, name, model, temper_time, temper_value);
                    end
                end
            end
            
            % to add the neighbour temperature
            for ii = 1:size(raw_neighbour,1)
                if isnan(raw_neighbour{ii,2})
                    %break
                else
                    number = raw_neighbour{ii,1};
                    name = raw_neighbour{ii,2};
                    model = raw_neighbour{ii,3};
                    temper_value = [];
                    if model == 0
                        temper_value = raw_neighbour{ii,4};
                        if sum(isnan(temper_value))
                            warning('Neighbour temperature with NaNs')
                            obj = obj.add_neighbour_set_temperature(building, number, name, model, ones(size(temper_value))*20);
                        else
                            obj = obj.add_neighbour_set_temperature(building, number, name, model, temper_value);
                        end
                    elseif model == 1
                        temper_value = [raw_neighbour{ii,4} raw_neighbour{ii,5}];
                        
                        if sum(isnan(temper_value))
                            warning('Neighbour temperature with NaNs')
                            obj = obj.add_neighbour_set_temperature(building, number, name, model, ones(size(temper_value))*20);
                        else
                            obj = obj.add_neighbour_set_temperature(building, number, name, model, temper_value);
                        end
                    elseif model == 2
                        temper_time = [0 31 59 90 120 151 181 212 243 273 304 334]*24*3600;
                        for ll = 1:12
                            temper_value(1,ll) = raw_neighbour{ii,ll+3};
                        end
                        
                        if sum(isnan(temper_value))
                            warning('Neighbour temperature with NaNs')
                            obj = obj.add_neighbour_set_temperature(building, number, name, model, temper_time, ones(size(temper_value))*20);
                        else
                            obj = obj.add_neighbour_set_temperature(building, number, name, model, temper_time, temper_value);
                        end
                    end
                end
            end
            
            % to add the weather
            name = raw_weather{1,2};
            name_txt = raw_weather{1,3};
            latitude = raw_weather{1,4};
            longitude = raw_weather{1,8};
            reference_meridian_for_time = raw_weather{1,12};
            obj = add_weather_from_file(obj, building, name_txt, name, latitude, longitude, reference_meridian_for_time);
        end
        
        function obj = ground_neighbour_weather_to_excel(obj, name_xls, building, variant_boundary)
            
            vollpfad = [pwd '\' name_xls];
            
            if exist(vollpfad, 'file')
                warning('Existing Excel file used.')
            else
                rootpath = fileparts(which('carnotUIBK'));
                copyfile([rootpath '\building_TEMPLATE.xlsx'], vollpfad);
                warning('New Excel file generated.')
            end
            
            % modify excel for write weather
            Excel = actxserver('Excel.Application');
            Excel.Workbooks.Open(vollpfad);
            warning('Excel file opened for writing! Do not interrupt this script!')
            
            % to delete the raws of the excel
            xlswrite1(name_xls,{''},'Boundary','B4:O9')
            xlswrite1(name_xls,{''},'Boundary','B14:O19')
            xlswrite1(name_xls,{''},'Boundary','B23:O23')
            
            count_ground = 0;
            count_neighbour = 0;
            
            for ii = 1:size(building.boundary(variant_boundary).ground,2)
                count_ground = count_ground+1;
                if strcmp(building.boundary(variant_boundary).ground(ii).name, 'none')
                    matrix_to_write_ground(count_ground,:) = [{''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}];
                else
                    ground_pfad = building.boundary(variant_boundary);
                    matrix_to_write_ground(count_ground,:) = [{ground_pfad.ground(count_ground).name} {2}...
                        {ground_pfad.ground(count_ground).temperature(1,2)} {ground_pfad.ground(count_ground).temperature(2,2)} {ground_pfad.ground(count_ground).temperature(3,2)}...
                        {ground_pfad.ground(count_ground).temperature(4,2)} {ground_pfad.ground(count_ground).temperature(5,2)} {ground_pfad.ground(count_ground).temperature(6,2)}...
                        {ground_pfad.ground(count_ground).temperature(7,2)} {ground_pfad.ground(count_ground).temperature(8,2)} {ground_pfad.ground(count_ground).temperature(9,2)}...
                        {ground_pfad.ground(count_ground).temperature(10,2)} {ground_pfad.ground(count_ground).temperature(11,2)} {ground_pfad.ground(count_ground).temperature(12,2)}];
                end
            end
            
            for ii = 1:size(building.boundary(variant_boundary).neighbour,2)
                count_neighbour = count_neighbour+1;
                 if strcmp(building.boundary(variant_boundary).neighbour(ii).name, 'none')
                    matrix_to_write_neigbour(count_neighbour,:) = [{''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}];
                 else
                    neighbour_pfad = building.boundary(variant_boundary);
                    matrix_to_write_neigbour(count_neighbour,:) = [{neighbour_pfad.neighbour(count_neighbour).name} {2}...
                        {neighbour_pfad.neighbour(count_neighbour).temperature(1,2)} {neighbour_pfad.neighbour(count_neighbour).temperature(2,2)} {neighbour_pfad.neighbour(count_neighbour).temperature(3,2)}...
                        {neighbour_pfad.neighbour(count_neighbour).temperature(4,2)} {neighbour_pfad.neighbour(count_neighbour).temperature(5,2)} {neighbour_pfad.neighbour(count_neighbour).temperature(6,2)}...
                        {neighbour_pfad.neighbour(count_neighbour).temperature(7,2)} {neighbour_pfad.neighbour(count_neighbour).temperature(8,2)} {neighbour_pfad.neighbour(count_neighbour).temperature(9,2)}...
                        {neighbour_pfad.neighbour(count_neighbour).temperature(10,2)} {neighbour_pfad.neighbour(count_neighbour).temperature(11,2)} {neighbour_pfad.neighbour(count_neighbour).temperature(12,2)}];
                 end
            end
            
            matrix_to_write_weather = [{building.boundary(variant_boundary).weather.name} {building.boundary(variant_boundary).weather.path} {building.boundary(variant_boundary).weather.latitude} {''} {''} {''} {building.boundary(variant_boundary).weather.longitude} {''} {''} {''} ...
                {building.boundary(variant_boundary).weather.latitude_timezone} {''} {''} {''}];
           
            xlswrite1(name_xls,matrix_to_write_ground,'Boundary',['B4:O' num2str(count_ground+3)]) 
            xlswrite1(name_xls,matrix_to_write_neigbour,'Boundary',['B14:O' num2str(count_neighbour+13)]) 
            xlswrite1(name_xls,matrix_to_write_weather,'Boundary','B23:O23') 
            
            
            Excel.ActiveWorkbook.Save
            Excel.Quit
            Excel.delete
            clear Excel
            warning('Excel file closed!')
            
        end
        
        function obj = add_weather_from_file(obj, building, name_txt, name, latitude, longitude, reference_meridian_for_time)
            weather = load(name_txt);
            
            timevec_weather = linspace(-8760*3600, (building.maxruntime+1)*8760*3600, ((building.maxruntime+2)*size(weather,1)-(building.maxruntime+1)))';
            weather_short = weather(1:(end-1),:);
            for jj = 1:(building.maxruntime+1)
                weather = [weather 
                                weather_short];
            end
            
            time_value = [timevec_weather timevec_weather]; 
            zenith = [timevec_weather weather(:,3) ];
            azimuth = [timevec_weather weather(:,4)];
            latitude_timezone = reference_meridian_for_time;
            radiation_beam_normal = [timevec_weather weather(:,5)];
            radiation_diffuse_horizontal = [timevec_weather weather(:,6)];
            t_ambient = [timevec_weather weather(:,7)];
            t_sky = [timevec_weather weather(:,8)];
            rh = [timevec_weather weather(:,9)];
            precip = [timevec_weather weather(:,10)];
            cloud = [timevec_weather weather(:,11)];
            p = [timevec_weather weather(:,12)];
            vw = [timevec_weather weather(:,13)];
            wdir = [timevec_weather weather(:,14)];
            incidence = [timevec_weather weather(:,15)];
            tetap = [timevec_weather weather(:,16)];
            tetas = [timevec_weather weather(:,17)];
            Idirect_surface = [timevec_weather weather(:,18)];
            Idiffuse_surface = [timevec_weather weather(:,19)];
            
            obj.weather = [obj.weather WEATHER(name, time_value, zenith, azimuth, latitude, longitude, latitude_timezone, radiation_beam_normal, radiation_diffuse_horizontal, t_ambient, t_sky, rh, precip, cloud, p, vw, wdir, incidence, tetap, tetas, Idirect_surface, Idiffuse_surface, name_txt)];
            
        end
        
        function obj = add_weather(obj, name, time_value, zenith, azimuth, latitude, longitude, latitude_timezone, radiation_beam_normal, radiation_diffuse_horizontal, t_ambient, t_sky, rh, precip, cloud, p, vw, wdir, incidence, tetap, tetas, Idirect_surface, Idiffuse_surface, path)
            obj.weather = [obj.weather WEATHER(name, time_value, zenith, azimuth, latitude, longitude, latitude_timezone, radiation_beam_normal, radiation_diffuse_horizontal, t_ambient, t_sky, rh, precip, cloud, p, vw, wdir, incidence, tetap, tetas, Idirect_surface, Idiffuse_surface, path)];
        end
        
        function obj = complete_ground_and_neighbour(obj, building)
            if length(obj.ground)<6
                len = length(obj.ground);
                name = 'none';
                model = obj.ground(len).model;
                temperature = obj.ground(len).temperature;
                for ii = 1 : (6-len)
                    obj.ground(len + ii) = GROUND();
                    obj.ground(len + ii).name = name;
                    obj.ground(len + ii).model = model;
                    obj.ground(len + ii).temperature = temperature;
                end
            end
            
            for ii = 1:length(obj.ground)
                if strcmp(obj.ground(ii).name, '')
                    obj = obj.add_ground_set_temperature(building, ii, 'none', 0, 20);
                end
            end
            
            if length(obj.neighbour)<6
                len = length(obj.neighbour);
                name = 'none';
                model = obj.neighbour(len).model;
                temperature= obj.neighbour(len).temperature;
                for ii = 1 : (6-len)
                    obj.neighbour(len + ii) = NEIGHBOUR();
                    obj.neighbour(len + ii).name = name;
                    obj.neighbour(len + ii).model = model;
                    obj.neighbour(len + ii).temperature = temperature;
                end
            end
            for ii = 1:length(obj.neighbour)
                if strcmp(obj.neighbour(ii).name, '')
                    obj = obj.add_neighbour_set_temperature(building, ii, 'none', 0, 20);
                end
            end
        end
        
        function neighbour = get_neighbour(obj, name)
        	ind = [];
            for jj = 1:length(obj.neighbour)
                if strcmp(obj.neighbour(jj).name,name)
                    ind = jj;
                    break
                end
            end
            %check if room is not existing
            if ind
                neighbour = obj.neighbour(ind);
            else
                error(['neighbour ' name ' not existing!'])
            end 
        end
        
        function weather = get_weather(obj, name)
        	ind = [];
            for jj = 1:length(obj.weather)
                if strcmp(obj.weather(jj).name,name)
                    ind = jj;
                    break
                end
            end
            %check if room is not existing
            if ind
                weather = obj.weather(ind);
            else
                error(['weather ' name ' not existing!'])
            end 
        end
        
        function ground = get_ground (obj, name)
        	ind = [];
            for jj = 1:length(obj.ground)
                if strcmp(obj.ground(jj).name,name)
                    ind = jj;
                    break
                end
            end
            %check if room is not existing
            if ind
                ground = obj.ground(ind);
            else
                error(['ground ' name ' not existing!'])
            end 
        end
        
    end
end

