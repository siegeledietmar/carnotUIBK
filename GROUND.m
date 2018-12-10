%% GROUND.m
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
% DS,EL     2017-03-12   initial revision v1.0

%%
classdef GROUND
    % GROUND
    
    properties
        
        name = '';
        model = '';                     % -1 no model
                                        % 0 fixed value for the whole year
                                        % 1 one fixed value for the
                                        % winter (first) and one fixed value
                                        % for the summer (second)
                                        % 2 sequence of values with time
                                        % (give it for one year)
        temperature = [];               % 2 columns (time and values)
        
    end
    
    methods
        function obj = GROUND(building, name, model, temperature_time, temperature_value)
            % to create a ground object
            % 1 ... building object
            % 2 ... name (for the user)
            % 3 ... model:    -1 no model
                            % 0 fixed value for the whole year
                            % 1 one fixed value for the
                            % winter (first) and one fixed value
                            % for the summer (second)
                            % 2 sequence of values with time
            % 4 ... temperature time: 1 or 2 value of temperature (model 0
            %                and 1) or time serie connected to value serie
            %               (temperature_value)
            % 5 ... optional series of temperature if temperature_time is a
            %               time serie
            if nargin==0
                
            elseif (nargin == 4)
                obj.name = name;
                obj.model = model;
                temperature_value = temperature_time;
                time_months = [0 31 59 90 120 151 181 212 243 273 304 334]*24*3600;
                temperature_time_ = [];
                for ii = -1:building.maxruntime
                    if ii == building.maxruntime
                        temperature_time_ = [temperature_time_ time_months+ii*365*24*3600 365*24*3600*(building.maxruntime+1)];
                    else
                        temperature_time_ = [temperature_time_ time_months+ii*365*24*3600];
                    end
                end
                if model == 0
                    temperature_value_ = temperature_value*ones(1,length(temperature_time_));
                end
                if model == 1
                    temperature_value_ = [];
                    for ii = -1:building.maxruntime
                        temperature_summer = temperature_value(1,2);
                        temperature_winter = temperature_value(1,1);
                        temperature_value_ = [temperature_value_ temperature_winter*ones(1,4) temperature_summer*ones(1,5) temperature_winter*ones(1,3)];
                    end
                    temperature_value_ = [temperature_value_ temperature_winter];
                end
                obj.temperature = [temperature_time_' temperature_value_'];
            elseif (nargin == 5)
                obj.name = name;
                obj.model = model;
                if model == 2
                    temperature_value_ = [];
                    temperature_time_ = [];
                    for ii = -1 : building.maxruntime
                        if ii == building.maxruntime
                            temperature_time_ = [temperature_time_ temperature_time+ii*365*24*3600 365*24*3600*(building.maxruntime+1)];
                            temperature_value_ = [temperature_value_ temperature_value temperature_value(1)];
                        else
                            temperature_time_ = [temperature_time_ temperature_time+ii*365*24*3600];
                            temperature_value_ = [temperature_value_ temperature_value];
                        end
                    end
                end
                obj.temperature = [temperature_time_' temperature_value_'];
            else
                error('Not a valid GROUND, too less or too much parameters!')
            end
        end
        
        function plot(obj)
            % to plot the ground temperature that are considered into
            % the building
            for ii = 1:length(obj)
                if strcmp(obj(ii).name, 'none')
                else
                    figure
                    name_plot = ['GROUND_' num2str(ii)];
                    gr = timeseries(obj(ii).temperature(:,2), obj(ii).temperature(:,1), 'Name', name_plot);
                    gr.TimeInfo.StartDate = '01-01-2014';
                    gr.TimeInfo.Format = 'dd-mmm';
                    gr.DataInfo.Units = '°C';
                    plot(gr, 'b-*')
                    grid on
                    xlabel('Time')
                    ylabel('Temperature (°C)')
                    title(name_plot)
                    xlim([datetime('31-Dec-2016 23:59:00') datetime('01-Jan-2018 00:00:00')])
                end
            end
        end
    end
end