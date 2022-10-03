%% GAINS.m
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

%%
classdef GAINS
    % GAINS
    
    properties
        name = 'new gains';
        gain = [];
        gain_to_zone = [];
    end
    
    methods
        function obj = GAINS()
            
        end
        
        function obj = add_gain(obj, name, model, type, timevalues, values1, values2, control)
            % to add a gain with the caracteristic
            % 1 ... name of the gain (for the user)
            % 2 ... model =  0...people, 1...light, 2...electricity, 3...moisture
            % 3 ... type =  % people: 0...W/ppl (const) 1...W/ppl (activity) 2...W/m²
                            % light: 0...W/light 1...W/m²
                            % electricity: 0...W/device 1...W/m²
                            % moisture: 0...kg/s
            % 4 ... time profile for the gains
            % 5 ... values1 =  for model 0, type 0 e 1: profile of people, for all the others: required profile
            % 6 ... values2 =  for model 0, type 0 e 1: activity profile
            % 7 ... control = 2; if gain is activated
            
            ind = [];
            for jj = 1:length(obj.gain)
                if strcmp(obj.gain(jj).name,name)
                    ind = strcmp(obj.gain(jj).name,name);
                    break
                end
            end
            if length(timevalues) == 2
                timevalues = [timevalues(1) (timevalues(end)-timevalues(1))/2 timevalues(end)];
                values1 = [values1 values1];
                values2 = [values2 values2];
            end
            
            % check if gain is already existing
            if ind
                error(['gain ' name ' already existing!'])
            else
                obj.gain = [obj.gain GAIN(name, model, type, timevalues, values1, values2, control)];
            end
        end
        
        function obj = assign_gain_to_zone(obj, name_gains, number_zone)
            % to assign the gain profile to the zone
            % 1 ... name of the gain profile (gain.name)
            % 2 ... number of the zone which you want to assign the profile to
            obj.gain_to_zone = [obj.gain_to_zone GAIN_TO_ZONE(name_gains, number_zone)];
        end
        
        function obj = gains_from_excel(obj, name_xls)
            % to take the gains profiles from the excel
            % 1 ... name of the excel file
            [~, ~, raw_gain] = xlsread(name_xls, 'Gains');
            
            raw_assigngaintozone = raw_gain(4:50,1:2);
            
            raw_person = raw_gain(4:33,4:31);
            raw_light = raw_gain(38:52,4:31);
            raw_electricity = raw_gain(57:71,4:31);
            raw_moisture = raw_gain(76:90,4:31);
            
            % to add the persons profiles
            for ii = 1:2:size(raw_person,1)
                if isnan(raw_person{ii,2})
                    break
                else
                    model = 0;
                    name = raw_person{ii,2};
                    type = raw_person{ii,3};
                    control = 2;
                    timevalues = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]*3600;
                    if type == 0 || type == 1
                        for ll = 1:24
                            values1(1,ll) = raw_person{ii,ll+4};
                            values2(1,ll) = raw_person{ii+1,ll+4};
                        end
                    end
                    if type == 2
                        for ll = 1:24
                            values1(1,ll) = raw_person{ii,ll+4};
                            values2(1,ll) = raw_person{ii,ll+4};
                        end
                    end
                    obj.gain = [obj.gain GAIN(name, model, type, timevalues, values1, values2, control)];
                end
            end
            
            % to add the light profiles
            for ii = 1:size(raw_light,1)
                if isnan(raw_light{ii,2})
                    break
                else
                    model = 1;
                    name = raw_light{ii,2};
                    type = raw_light{ii,3};
                    control = 2;
                    timevalues = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]*3600;
                    if type == 0 || type == 1
                        for ll = 1:24
                            values1(1,ll) = raw_light{ii,ll+4};
                            values2(1,ll) = 0;
                        end
                    end
                    obj.gain = [obj.gain GAIN(name, model, type, timevalues, values1, values2, control)];
                end
            end
            
            % to add the electricity profiles
            for ii = 1:size(raw_electricity,1)
                if isnan(raw_electricity{ii,2})
                    break
                else
                    model = 2;
                    name = raw_electricity{ii,2};
                    type = raw_electricity{ii,3};
                    control = 2;
                    timevalues = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]*3600;
                    if type == 0 || type == 1
                        for ll = 1:24
                            values1(1,ll) = raw_electricity{ii,ll+4};
                            values2(1,ll) = 0;
                        end
                    end
                    obj.gain = [obj.gain GAIN(name, model, type, timevalues, values1, values2, control)];
                end
            end
            
            % to add the moisture profiles
            for ii = 1:size(raw_moisture,1)
                if isnan(raw_moisture{ii,2})
                    break
                else
                    model = 3;
                    name = raw_moisture{ii,2};
                    type = raw_moisture{ii,3};
                    control = 2;
                    timevalues = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]*3600;
                    if type == 0 || type == 1
                        for ll = 1:24
                            values1(1,ll) = raw_moisture{ii,ll+4};
                            values2(1,ll) = 0;
                        end
                    end
                    obj.gain = [obj.gain GAIN(name, model, type, timevalues, values1, values2, control)];
                end
            end
            
            % to assign the gain profile to the zone
            for ii = 1:size(raw_assigngaintozone,1)
                if isnan(raw_assigngaintozone{ii,1})
                    break
                else
                    name_gains = raw_assigngaintozone{ii,1};
                    number_zone = raw_assigngaintozone{ii,2};
                end
                obj.gain_to_zone = [obj.gain_to_zone GAIN_TO_ZONE(name_gains, number_zone)];
            end
        end
        
        function obj = gains_to_excel(obj, name_xls, building, variant_gains)
            
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
            xlswrite1(name_xls,{''},'Gains','A4:B100')
            
            xlswrite1(name_xls,{''},'Gains','E4:F33')
            xlswrite1(name_xls,{''},'Gains','H4:AE33')
            
            xlswrite1(name_xls,{''},'Gains','E38:F52')
            xlswrite1(name_xls,{''},'Gains','H38:AE52')
            
            xlswrite1(name_xls,{''},'Gains','E57:F71')
            xlswrite1(name_xls,{''},'Gains','H57:AE71')
            
            xlswrite1(name_xls,{''},'Gains','E76:F90')
            xlswrite1(name_xls,{''},'Gains','H76:AE90')
            
            
            count_people = 0;
            count_light = 0;
            count_electr = 0;
            count_moist = 0;
            
            for ii = 1:size(building.gains(variant_gains).gain,2)
                model = building.gains(variant_gains).gain(ii).model;
                
                switch model
                    case 0 % people
                        count_people = count_people+1;
                        
                        name_people{count_people*2-1} = building.gains(variant_gains).gain(ii).name;
                        type_people{count_people*2-1} = building.gains(variant_gains).gain(ii).type;

                        if type_people{count_people*2-1} == 2
                            if building.gains(variant_gains).gain(ii).values1 == 24
                                values1_people(count_people*2-1,1:24) = building.gains(variant_gains).gain(ii).values1;
                            else
                                values1_people(count_people*2-1,1:24) = building.gains(variant_gains).gain(ii).values1(1);
                            end
                            matrix_to_write_people(count_people*2-1,:) = [{building.gains(variant_gains).gain(ii).name} {building.gains(variant_gains).gain(ii).type} {'values 1'} {values1_people(count_people*2-1,1)} {values1_people(count_people*2-1,2)} {values1_people(count_people*2-1,3)} {values1_people(count_people*2-1,4)} {values1_people(count_people*2-1,5)} {values1_people(count_people*2-1,6)} {values1_people(count_people*2-1,7)} {values1_people(count_people*2-1,8)} {values1_people(count_people*2-1,9)} {values1_people(count_people*2-1,10)} {values1_people(count_people*2-1,11)} {values1_people(count_people*2-1,12)} {values1_people(count_people*2-1,13)} {values1_people(count_people*2-1,14)} {values1_people(count_people*2-1,15)} {values1_people(count_people*2-1,16)} {values1_people(count_people*2-1,17)} {values1_people(count_people*2-1,18)} {values1_people(count_people*2-1,19)} {values1_people(count_people*2-1,20)} {values1_people(count_people*2-1,21)} {values1_people(count_people*2-1,22)} {values1_people(count_people*2-1,23)} {values1_people(count_people*2-1,24)}];
                            matrix_to_write_people(count_people*2,:) = [{''} {''} {'values 2'} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}];
                        else
                            if building.gains(variant_gains).gain(ii).values1 == 24
                                values1_people(count_people*2-1,1:24) = building.gains(variant_gains).gain(ii).values1(1:24);
                                values1_people(count_people*2,1:24) = building.gains(variant_gains).gain(ii).values2(1:24);
                            else
                                values1_people(count_people*2-1,1:24) = building.gains(variant_gains).gain(ii).values1(1);
                                values1_people(count_people*2,1:24) = building.gains(variant_gains).gain(ii).values2(1);
                            end
                            
                            matrix_to_write_people(count_people*2-1,:) = [{building.gains(variant_gains).gain(ii).name} {building.gains(variant_gains).gain(ii).type} {'values 1'} {values1_people(count_people*2-1,1)} {values1_people(count_people*2-1,2)} {values1_people(count_people*2-1,3)} {values1_people(count_people*2-1,4)} {values1_people(count_people*2-1,5)} {values1_people(count_people*2-1,6)} {values1_people(count_people*2-1,7)} {values1_people(count_people*2-1,8)} {values1_people(count_people*2-1,9)} {values1_people(count_people*2-1,10)} {values1_people(count_people*2-1,11)} {values1_people(count_people*2-1,12)} {values1_people(count_people*2-1,13)} {values1_people(count_people*2-1,14)} {values1_people(count_people*2-1,15)} {values1_people(count_people*2-1,16)} {values1_people(count_people*2-1,17)} {values1_people(count_people*2-1,18)} {values1_people(count_people*2-1,19)} {values1_people(count_people*2-1,20)} {values1_people(count_people*2-1,21)} {values1_people(count_people*2-1,22)} {values1_people(count_people*2-1,23)} {values1_people(count_people*2-1,24)}];
                            matrix_to_write_people(count_people*2,:) = [{''} {''} {'values 2'} {values1_people(count_people*2,1)} {values1_people(count_people*2,2)} {values1_people(count_people*2,3)} {values1_people(count_people*2,4)} {values1_people(count_people*2,5)} {values1_people(count_people*2,6)} {values1_people(count_people*2,7)} {values1_people(count_people*2,8)} {values1_people(count_people*2,9)} {values1_people(count_people*2,10)} {values1_people(count_people*2,11)} {values1_people(count_people*2,12)} {values1_people(count_people*2,13)} {values1_people(count_people*2,14)} {values1_people(count_people*2,15)} {values1_people(count_people*2,16)} {values1_people(count_people*2,17)} {values1_people(count_people*2,18)} {values1_people(count_people*2,19)} {values1_people(count_people*2,20)} {values1_people(count_people*2,21)} {values1_people(count_people*2,22)} {values1_people(count_people*2,23)} {values1_people(count_people*2,24)}];
                        end
                        
                    case 1 % light
                        count_light = count_light+1;
                        values1_light(count_light,1:24) = building.gains(variant_gains).gain(ii).values1;
                        matrix_to_write_light(count_light,:) = [{building.gains(variant_gains).gain(ii).name} {building.gains(variant_gains).gain(ii).type} {'values 1'} {values1_light(count_light,1)} {values1_light(count_light,2)} {values1_light(count_light,3)} {values1_light(count_light,4)} {values1_light(count_light,5)} {values1_light(count_light,6)} {values1_light(count_light,7)} {values1_light(count_light,8)} {values1_light(count_light,9)} {values1_light(count_light,10)} {values1_light(count_light,11)} {values1_light(count_light,12)} {values1_light(count_light,13)} {values1_light(count_light,14)} {values1_light(count_light,15)} {values1_light(count_light,16)} {values1_light(count_light,17)} {values1_light(count_light,18)} {values1_light(count_light,19)} {values1_light(count_light,20)} {values1_light(count_light,21)} {values1_light(count_light,22)} {values1_light(count_light,23)} {values1_light(count_light,24)}];
                        
                    case 2 % electricity
                        count_electr = count_electr+1;
                        values1_electr(count_electr,1:24) = building.gains(variant_gains).gain(ii).values1;
                        matrix_to_write_electr(count_electr,:) = [{building.gains(variant_gains).gain(ii).name} {building.gains(variant_gains).gain(ii).type} {'values 1'} {values1_electr(count_electr,1)} {values1_electr(count_electr,2)} {values1_electr(count_electr,3)} {values1_electr(count_electr,4)} {values1_electr(count_electr,5)} {values1_electr(count_electr,6)} {values1_electr(count_electr,7)} {values1_electr(count_electr,8)} {values1_electr(count_electr,9)} {values1_electr(count_electr,10)} {values1_electr(count_electr,11)} {values1_electr(count_electr,12)} {values1_electr(count_electr,13)} {values1_electr(count_electr,14)} {values1_electr(count_electr,15)} {values1_electr(count_electr,16)} {values1_electr(count_electr,17)} {values1_electr(count_electr,18)} {values1_electr(count_electr,19)} {values1_electr(count_electr,20)} {values1_electr(count_electr,21)} {values1_electr(count_electr,22)} {values1_electr(count_electr,23)} {values1_electr(count_electr,24)}];
                        
                    case 3 % moisture
                        count_moist = count_moist+1;
                        values1_moist(count_moist,1:24) = building.gains(variant_gains).gain(ii).values1;
                        matrix_to_write_moist(count_moist,:) = [{building.gains(variant_gains).gain(ii).name} {building.gains(variant_gains).gain(ii).type} {'values 1'} {values1_moist(count_moist,1)} {values1_moist(count_moist,2)} {values1_moist(count_moist,3)} {values1_moist(count_moist,4)} {values1_moist(count_moist,5)} {values1_moist(count_moist,6)} {values1_moist(count_moist,7)} {values1_moist(count_moist,8)} {values1_moist(count_moist,9)} {values1_moist(count_moist,10)} {values1_moist(count_moist,11)} {values1_moist(count_moist,12)} {values1_moist(count_moist,13)} {values1_moist(count_moist,14)} {values1_moist(count_moist,15)} {values1_moist(count_moist,16)} {values1_moist(count_moist,17)} {values1_moist(count_moist,18)} {values1_moist(count_moist,19)} {values1_moist(count_moist,20)} {values1_moist(count_moist,21)} {values1_moist(count_moist,22)} {values1_moist(count_moist,23)} {values1_moist(count_moist,24)}];
                end                
            end
            if exist('matrix_to_write_people') ~=  0
               xlswrite1(name_xls,matrix_to_write_people,'Gains',['E4:AE' num2str(count_people*2+3)]) 
            end
            if exist('matrix_to_write_light') ~=  0
               xlswrite1(name_xls,matrix_to_write_light,'Gains',['E38:AE' num2str(count_light+37)])
            end
            if exist('matrix_to_write_electr') ~=  0
               xlswrite1(name_xls,matrix_to_write_electr,'Gains',['E57:AE' num2str(count_electr+56)]) 
            end
            if exist('matrix_to_write_moist') ~=  0
               xlswrite1(name_xls,matrix_to_write_moist,'Gains',['E76:AE' num2str(count_moist+75)])
            end
            
            count_gaintozone = 0;
            
            for jj=1:size(building.gains(variant_gains).gain_to_zone,2)
                count_gaintozone = count_gaintozone +1;
                matrix_to_write_gaintozone(count_gaintozone,:) = [{building.gains(variant_gains).gain_to_zone(jj).name_gains} {building.gains(variant_gains).gain_to_zone(jj).number_of_zone}];
            end
            
            xlswrite1(name_xls,matrix_to_write_gaintozone,'Gains',['A4:B' num2str(count_gaintozone+3)])
            
            Excel.ActiveWorkbook.Save
            Excel.Quit
            Excel.delete
            clear Excel
            warning('Excel file closed!')
        end
        
        function gain = get_gain(obj, name)
            % to get the gain
            % 1 ... name of the gain
            ind = [];
            for jj = 1:length(obj.gain)
                if strcmp(obj.gain(jj).name,name)
                    ind = jj;
                    break
                end
            end
            
            if ind
                gain = obj.gain(ind);
            else
                error(['gain ' name ' not existing!'])
            end
        end
        
        function obj = gain_from_PHPP(obj, name_xls_PHPP, language, version)
            % to take the gains from the PHPP
            % 1 ... name of the file xls of the PHPP
            % 2 ... language: 1:german, 0:english
            % 3 ... version of PHPP used
            filename = name_xls_PHPP;
            
            switch version
                case '9.1'
                    if language
                        sheet = 'Nachweis';
                    else
                        sheet = 'Verification';
                    end
                    range = 'I28:K34';
                    [~, ~, data] = xlsread(filename, sheet, range);
                    index_raw_intloads = 1;
                    index_column_intloads = 3;
                    index_raw_area = 7;
                    index_column_area = 1;
                    
               case '10.2'
                if language
                    sheet = 'Nachweis';
                else
                    sheet = 'Verification';
                end
                range = 'I29:K35';
                [~, ~, data] = xlsread(filename, sheet, range);
                index_raw_intloads = 1;
                index_column_intloads = 3;
                index_raw_area = 7;
                index_column_area = 1;
            end
            timevalues = [0:24]*3600;
            values1 = ones(1,24)*data{index_raw_intloads,index_column_intloads};
            values2 = values1; % not used
            control = 2;
            obj = obj.add_gain('Person_W/m²_from_PHPP', 0, 2, timevalues, values1, values2, control);
        end
        
        function plot_gain(obj)
            % to plot all the gain profiles
            for ii = 1:length(obj.gain)
                figure
                %stairs(obj.gain(ii).timevalues/3600,[obj.gain(ii).values1 obj.gain(ii).values1(1)])
                xxx = [obj.gain(ii).timevalues/3600;obj.gain(ii).timevalues/3600];
                yyy = [obj.gain(ii).values1 obj.gain(ii).values1(1); obj.gain(ii).values1 obj.gain(ii).values1(1)];
                h = area(xxx([2:end end]),yyy(1:end));
                h.FaceColor = [0 0.75 0.75];
                name_mod = strrep(obj.gain(ii).name,'_',' ');
                title(name_mod)
                xlim([0 24])
                set(gca,'XTick',1:24)
                if obj.gain(ii).model ==0 && (obj.gain(ii).type == 0 || obj.gain(ii).type == 1)
                    ylabel('People (people)')
                elseif obj.gain(ii).model == 0 && obj.gain(ii).type == 2
                    ylabel('People (W/m²)')
                elseif obj.gain(ii).model == 1 && obj.gain(ii).type == 0
                    ylabel('Light (W/light)')
                elseif obj.gain(ii).model == 1 && obj.gain(ii).type == 1
                    ylabel('Light (W/m²)')
                elseif obj.gain(ii).model == 2 && obj.gain(ii).type == 0
                    ylabel('Electricity (W/device)')
                elseif obj.gain(ii).model == 2 && obj.gain(ii).type == 1
                    ylabel('Electricity (W/m²)')
                elseif obj.gain(ii).model == 3 && obj.gain(ii).type == 0
                    ylabel('Moisture (kg/s)')
                end
                xlabel('Time (h)')
            end
        end
        
        function calculate_gain(obj, build, thzone)
            % to plot all the gain profiles
            energy = 0;
            area_tot = 0;
            zone_intgains_ = [];
            zone_intgains(1:10) = 0;
            for ii = 1:length(obj.gain)
                energy_year = 0;
                area_zone = [];
                [energy_year area_zone zone_intgains] = calc_gains_ppl(obj, build, ii, obj.gain(ii).model, obj.gain(ii).type, zone_intgains, thzone);
                energy = energy+sum(energy_year);
            end
            for mm = 1:length(build.thermalzone(thzone).zone)
                if zone_intgains(mm) == 1
                    area_tot = area_tot + build.thermalzone(thzone).zone(mm).heated_area;
                end
            end
            disp(' ')
            disp(['Energy year tot: ' num2str(energy) ' kWh/y'])
            disp(['Specific Energy: ' num2str(energy*1000/(area_tot*24*365)) ' W/m²'])
        end
    end
    
    methods (Access = private)
        function [energy_year area_zone zone_intgains energy_day  energy_specific] = calc_gains_ppl(obj, build, n_gain, model, type, zone_intgains, thzone)
            if obj.gain(n_gain).model ==0 && (obj.gain(n_gain).type == 0 || obj.gain(n_gain).type == 1)
                disp(' ')
                disp('People (people)')
            elseif obj.gain(n_gain).model == 0 && obj.gain(n_gain).type == 2
                disp(' ')
                disp('People (W/m²)')
            elseif obj.gain(n_gain).model == 1 && obj.gain(n_gain).type == 0
                disp(' ')
                disp('Light (W/light)')
            elseif obj.gain(n_gain).model == 1 && obj.gain(n_gain).type == 1
                disp(' ')
                disp('Light (W/m²)')
            elseif obj.gain(n_gain).model == 2 && obj.gain(n_gain).type == 0
                disp(' ')
                disp('Electricity (W/device)')
            elseif obj.gain(n_gain).model == 2 && obj.gain(n_gain).type == 1
                disp(' ')
                disp('Electricity (W/m²)')
            elseif obj.gain(n_gain).model == 3 && obj.gain(n_gain).type == 0
                disp(' ')
                disp('Moisture (kg/s)')
            end
            zone_number = [];
            count_zone = 0;
            area_zone = [];
            for iii = 1:length(obj.gain_to_zone)
                if strcmp(obj.gain(n_gain).name, obj.gain_to_zone(iii).name_gains)
                    count_zone = count_zone + 1;
                    zone_intgains(obj.gain_to_zone(iii).number_of_zone) = 1;
                    zone_number(count_zone) = obj.gain_to_zone(iii).number_of_zone;
                    area_zone(count_zone) = build.thermalzone(thzone).zone(obj.gain_to_zone(iii).number_of_zone).heated_area;
                end
            end
            if (model == 0 && (type ==0 || type ==1)) || (model == 1 && type ==0) || (model == 2 && type ==0) || model == 3
                energy_day = 0;
                energy_year = 0;
                if model == 0 && (type ==0 || type ==1)
                    watt = 60;
                elseif (model == 1 && type ==0) || (model == 2 && type ==0) 
                    watt = 1;
                elseif model == 3
                    watt = -2500*1000;
                end
                for jj = 1:24
                    energy_day = energy_day + watt*(obj.gain(n_gain).values1(jj)*obj.gain(n_gain).timeduration(jj)); % W*s/ day
                end
                energy_year = energy_day*365/(1000*3600); %kWh/y
                disp([obj.gain(n_gain).name ', zone(s): ' num2str(zone_number(1:end))])
                disp(['Energy year = : ' num2str(energy_year) ' kWh/y'])
                specific_energy = [];
                for kk = 1:length(zone_number)
                    specific_energy(kk) = energy_day/(24*3600*area_zone(kk)); %W/m²
                    disp(['Specific Energy, zone ' num2str(zone_number(kk)) ' = ' num2str(specific_energy(kk)) ' W/m²'])
                end
            else
                disp([obj.gain(n_gain).name ', zone(s): ' num2str(zone_number(1:end))])
                specific_energy = [];
                energy_year = [];
                for kk = 1:length(zone_number)
                    energy_day = 0;
                    for jj = 1:24
                        energy_day = energy_day + (obj.gain(n_gain).values1(jj)*obj.gain(n_gain).timeduration(jj))*area_zone(kk); % W*s/ day
                    end
                    specific_energy(kk) = energy_day/(24*3600*area_zone(kk)); %W/m²
                    energy_year(kk) = energy_day*365/(1000*3600); %kWh/y
                    disp(['Energy year, zone ' num2str(zone_number(kk)) ' = ' num2str(energy_year(kk)) ' kWh/y'])
                    disp(['Specific Energy, zone ' num2str(zone_number(kk)) ' = ' num2str(specific_energy(kk)) ' W/m²'])
                end
            end     
        end
    end
end

