%% STRUCTURE.m
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
% EL        2018-05-25   v1.3: prepare for the import from Excel
% DS,EL	    2019-01-24   updates for GUI v2.0

%%
classdef STRUCTURE
    % STRUCTURE
    % heat flux always from 1 to 2!
    % 1 ... inside
    % 2 ... outside
    
    properties
        name = '';                  % name of structure
        type = 0;                   % type of the structure: 0 ... wall, 1 ... window, 2 ... thermal bridge
        model = -1;                 % wall: 0 ... UA, 1 ... RC, 2 ... hygrothermal
                                    % window: 0 ... TF, 1 ... CARNOT, 2 ... TF shading
                                    % thermal bridge: 0 ... UA, 1 ... RC
        parameter = [];             % parameter set
        emission_1 = 0.94;          % emission inside / [-]
        emission_2 = 0.60;          % emission outside / [-]
        absorption_2 = 0.65;        % absorption outside / [-]
    end
    
    methods
        function obj = STRUCTURE(name, parameter, emission_1, emission_2, absorption_2)
            if (nargin < 2) || (nargin > 5)
                error('Not a valid STRUCTURE!')
            elseif nargin < 3
                obj.name = name;
                obj.parameter = load(parameter);
                obj.type = obj.parameter.type;
                obj.model = obj.parameter.model;
            else
                obj.name = name;
                if isstruct(parameter) %30_01_2018
                    obj.parameter = parameter;
                else
                    if ~isempty(parameter)
                        try
                            obj.parameter = load(parameter);
                        catch
                            try
                                obj.parameter = parameter;
                            catch
                                obj.parameter = [];
                            end
                        end
                    else
                        obj.parameter = [];
                    end
                end
                obj.type = obj.parameter.type;
                obj.model = obj.parameter.model;
                obj.emission_1 = emission_1;
                obj.emission_2 = emission_2;
                obj.absorption_2 = absorption_2;
            end
        end
        
        function obj = set_parameter(obj,parameter)
            % to load the file of the parameters and set the parameters of
            % the structure
            % 1 ... name of the parameter file.mat
            obj.parameter = load(parameter);
            obj.model = obj.parameter.model;
        end
        
        function plot(obj, number_str)
            % to plot the structures
            % number of structure that you want to plot
            if nargin == 1
                for ii = 1:length(obj)
                    % check if d is existing
                    if isfield(obj(ii).parameter,'d')  
                        if length(obj(ii).parameter.d) == 1 && (obj(ii).parameter.d) == 1
                        else
                            shape = zeros(length(obj(ii).parameter.d)+1,1);
                            for jj = 1:length(obj(ii).parameter.d)
                                shape(jj+1) = shape(jj) + obj(ii).parameter.d(jj);
                            end
                            % open figure
                                f1 = figure; hold all
                                % if colors are existing
                                if isfield(obj(ii).parameter,'layers_colors')
                                    for jj = 1:length(obj(ii).parameter.d)
                                        rectangle('Position',[shape(jj),-0.1,obj(ii).parameter.d(jj),1.2],'FaceColor',obj(ii).parameter.layers_colors(jj,:),'EdgeColor','k','LineWidth',2)
                                    end
                                else
                                    for jj = 1:length(obj(ii).parameter.d)
                                        rectangle('Position',[shape(jj),-0.1,obj(ii).parameter.d(jj),1.2],'FaceColor',rand(1,3),'EdgeColor','k','LineWidth',2)
                                    end
                                end

                                % beuken model
                                if obj(ii).model>=1 && isfield(obj(ii).parameter,'xmesh_beu')
                                    red_point = 0;
                                    for jj=1:length(obj(ii).parameter.xmesh_beu)
                                        red_point = red_point + obj(ii).parameter.xmesh_beu(jj);
                                        plot(red_point,0.5,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor',[0.1,0.1,0.1],'MarkerSize',10)
                                    end
                                    plot(0,0.5,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor',[0.1,0.1,0.1],'MarkerSize',10)
                                    plot(shape(end),0.5,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor',[0.1,0.1,0.1],'MarkerSize',10)  
                                    plot([0:0.1:shape(end) shape(end)],ones(length(0:0.1:shape(end))+1,1)*0.5,'color','r','LineStyle','-')
                                end

                                % hygrothermal model
                                if obj(ii).model==2 && isfield(obj(ii).parameter,'xmesh_hygro')
                                    for jj=1:length(obj(ii).parameter.xmesh_hygro)
                                        plot(ones(1,2)*obj(ii).parameter.xmesh_hygro(jj),[-0.1 1.1],'color',[0.5,0.5,0.5],'LineStyle','--')
                                    end
                                end

                                % if names exist
                                if isfield(obj(ii).parameter,'layers_names')
                                    for jj=1:length(obj(ii).parameter.layers_names)
                                        if obj(ii).parameter.d(jj) < 0.01
                                            text(shape(jj)-obj(ii).parameter.d(jj),0.2,obj(ii).parameter.layers_names(jj),'HorizontalAlignment','center','Rotation',90)
                                        else
                                            text(shape(jj)+obj(ii).parameter.d(jj)/2,0.2,obj(ii).parameter.layers_names(jj),'HorizontalAlignment','center','Rotation',90)
                                        end
                                    end
                                end

                                xlim([-.05 shape(end)+0.05])
                                ylim([0 1])
                                name_mod = strrep(obj(ii).name,'_',' ');

                                title(name_mod)
                        end
                    else
                        warning('No dimensions available to plot!')
                    end
                end
            elseif nargin == 2
                % check if d is existing
                if isfield(obj(number_str).parameter,'d')
                    if length(obj(ii).parameter.d) == 1 && (obj(ii).parameter.d) == 1
                    else
                        shape = zeros(length(obj(number_str).parameter.d)+1,1);
                        for jj = 1:length(obj(number_str).parameter.d)
                            shape(jj+1) = shape(jj) + obj(number_str).parameter.d(jj);
                        end

                        % open figure
                        figure, hold all
                            % if colors are existing
                            if isfield(obj(number_str).parameter,'layers_colors')
                                for jj = 1:length(obj(number_str).parameter.d)
                                    rectangle('Position',[shape(jj),-0.1,obj(number_str).parameter.d(jj),1.2],'FaceColor',obj(number_str).parameter.layers_colors(jj,:),'EdgeColor','k','LineWidth',2)
                                end
                            else
                                for jj = 1:length(obj(number_str).parameter.d)
                                    rectangle('Position',[shape(jj),-0.1,obj(number_str).parameter.d(jj),1.2],'FaceColor',rand(1,3),'EdgeColor','k','LineWidth',2)
                                end
                            end

                            % beuken model
                            if obj(number_str).model==1 && isfield(obj(number_str).parameter,'xmesh_beu')
                                red_point = 0;
                                for jj=1:length(obj(number_str).parameter.xmesh_beu)
                                    red_point = red_point + obj(number_str).parameter.xmesh_beu(jj);
                                    plot(red_point,0.5,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor',[0.1,0.1,0.1],'MarkerSize',10) 
                                    obj(number_str).parameter.xmesh_beu(jj)
                                end
                                plot(0,0.5,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor',[0.1,0.1,0.1],'MarkerSize',10)
                                plot(shape(end),0.5,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor',[0.1,0.1,0.1],'MarkerSize',10)  
                                plot([0:0.1:shape(end) shape(end)],ones(length(0:0.1:shape(end))+1,1)*0.5,'color','r','LineStyle','-')
                            end

                            % hygrothermal model
                            if obj(number_str).model==2 && isfield(obj(number_str).parameter,'xmesh_hygro')
                                for jj=1:length(obj(number_str).parameter.xmesh_hygro)
                                    plot(ones(1,2)*obj(number_str).parameter.xmesh_hygro(jj),[-0.1 1.1],'color',[0.5,0.5,0.5],'LineStyle','--')
                                end
                            end

                            % if names exist
                            if isfield(obj(number_str).parameter,'layers_names')
                                for jj=1:length(obj(number_str).parameter.layers_names)
                                    if obj(number_str).parameter.d(jj) < 0.01
                                        text(shape(jj)-obj(number_str).parameter.d(jj),0.2,obj(number_str).parameter.layers_names(jj),'HorizontalAlignment','center','Rotation',90)
                                    else
                                        text(shape(jj)+obj(number_str).parameter.d(jj)/2,0.2,obj(number_str).parameter.layers_names(jj),'HorizontalAlignment','center','Rotation',90)
                                    end
                                end
                            end

                            xlim([-.05 shape(end)+0.05])
                            ylim([0 1])
                            name_mod = strrep(obj(number_str).name,'_',' ');

                            title(name_mod)
                    end
                else
                    warning('No dimensions available to plot!')
                end
            end
        end
        
        function disp(obj, disp_more)
            % to display all the structure with the model and the U value,
            % 1 ... optional: only if you want to display also other
            % parameters (emission, absorption, layers)
            for ii = 1:length(obj)
                disp([obj(ii).name])
                if obj(ii).type == 0
                    if obj(ii).model == 0
                        disp('Wall model: UA')
                    elseif obj(ii).model == 1
                        disp('Wall model: RC')
                    elseif obj(ii).model == 2
                        disp('Wall model: Hygrothermal')
                    end
                elseif obj(ii).type == 1
                    if obj(ii).model == 0
                        disp('Window model: TF')
                    elseif obj(ii).model == 1
                        disp('Window model: CARNOT')
                    elseif obj(ii).model == 2
                        disp('Window model: TF shading')
                    end
                elseif obj(ii).type == 2
                    if obj(ii).model == 0
                        disp('Thermal bridge model: UA')
                    elseif obj(ii).model == 1
                        disp('Thermal bridge model: RC')
                    end
                end
                
                if obj(ii).type == 0 || obj(ii).type == 2
                    disp(['L value (U value without heat transfer coefficients): ' num2str(obj(ii).parameter.U) ' W/(m².K)'])
                else
                    disp(['U value glass: ' num2str(obj(ii).parameter.U_g) ' W/(m²K),  U value frame: ' num2str(obj(ii).parameter.U_f) ' W/(m²K)'])
                end
%                 if nargin == 2
                    if obj(ii).type == 0 || obj(ii).type == 2
                        disp(['Emission: ' num2str(obj(ii).emission_1) ', Emission: ' num2str(obj(ii).emission_2) ', Absorbtion: ' num2str(obj(ii).absorption_2) ])
                        if isfield(obj(ii).parameter, 'tau')
                            disp(['Tau: ' num2str(obj(ii).parameter.tau) ])
                        end
                    else
                        disp(['Emission: ' num2str(obj(ii).emission_1) ', Emission: ' num2str(obj(ii).emission_2) ', Absorbtion: ' num2str(obj(ii).absorption_2) ', Tau: ' num2str(obj(ii).parameter.tau_g_w) ])
                    end
                    if obj(ii).type == 0 || obj(ii).type == 2
                        if length(obj(ii).parameter.d) ==1 && obj(ii).parameter.d==1
                        else
                            disp(['Layers: '])
                            for ll = 1:length(obj(ii).parameter.d)
                                disp([obj(ii).parameter.layers_names{ll} ' -> ' num2str(obj(ii).parameter.d(ll)*100) ' cm']);
                            end
                        end
                    end
%                 end
                disp('   ')
            end
        end
        
        function disp_detailed(obj, num_str, disp_more)
            % to display one structure with the model and the U value,
            % 1 ... number of the structure that you want to display
            % 2 ... optional: only if you want to display also other
            % parameters (emission, absorption, layers)
            disp([obj(num_str).name])
            if obj(num_str).type == 0
                if obj(num_str).model == 0
                    disp('Wall model: UA')
                elseif obj(num_str).model == 1
                    disp('Wall model: RC')
                elseif obj(num_str).model == 2
                    disp('Wall model: Hygrothermal')
                end
            elseif obj(num_str).type == 1
                if obj(num_str).model == 0
                    disp('Window model: TF')
                elseif obj(num_str).model == 1
                    disp('Window model: CARNOT')
                elseif obj(num_str).model == 2
                    disp('Window model: TF shading')
                end
            elseif obj(num_str).type == 2
                if obj(num_str).model == 0
                    disp('Thermal bridge model: UA')
                elseif obj(num_str).model == 1
                    disp('Thermal bridge model: RC')
                end
            end

            if obj(num_str).type == 0 || obj(num_str).type == 2
                disp(['L value (U value without heat transfer coefficients): ' num2str(obj(num_str).parameter.U) ' W/(m²K)'])
            else
                disp(['U value glass: ' num2str(obj(num_str).parameter.U_g) ' W/(m²K),  U value frame: ' num2str(obj(num_str).parameter.U_f) ' W/(m²K)'])
            end
            if nargin == 3
                if obj(num_str).type == 0 || obj(num_str).type == 2
                    disp(['Emission: ' num2str(obj(num_str).emission_1) ', Emission: ' num2str(obj(num_str).emission_2) ', Absorbtion: ' num2str(obj(num_str).absorption_2) ])
                    if isfield(obj(num_str).parameter, 'tau')
                        disp(['Tau: ' num2str(obj(num_str).parameter.tau) ])
                    end
                else
                    disp(['Emission: ' num2str(obj(num_str).emission_1) ', Emission: ' num2str(obj(num_str).emission_2) ', Absorbtion: ' num2str(obj(num_str).absorption_2) ', Tau: ' num2str(obj(num_str).parameter.tau_g_w) ])
                end
                if obj(num_str).type == 0 || obj(num_str).type == 2
                    if length(obj(num_str).parameter.d) ==1 && obj(num_str).parameter.d==1
                    else
                        disp(['Layers: '])
                        for ll = 1:length(obj(num_str).parameter.d)
                            disp([obj(num_str).parameter.layers_names{ll} ' -> ' num2str(obj(num_str).parameter.d(ll)*100) ' cm']);
                        end
                    end
                end
            end
            disp('   ') 
        end
        
    end
end


