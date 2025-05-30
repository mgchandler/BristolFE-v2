function h = fn_display_result_v2(nds, els, display_options)
%SUMMARY
%   Displays mesh from 2D model, returning handle to patches for later
%   animations
%USAGE
%   fn_display_result(nodes, elements, display_options) to display mesh OR
%INPUTS
%   nodes - n x 2 matrix of nodal coordinates. The row number is the node
%   number; columns 1 and 2 are the x and y coordinates of the node.
%   elements - m x 3 matrix of element nodes. The row number is the element
%   number; columns 1, 2 and 3 are the node numbers of the 3 nodes for each
%   triangular element
%   display_options - structured variable allowing optional plotting properties to
%   be set. See below for defaults. In particular:
%   default_options.node_sets_to_plot - allows specific nodes to be plotted
%   in a particular color. It is a vector of structured variables with
%   fields nd and col. nd is a vector of node indices and col is the color
%   (e.g. 'r') in which nodes in that set will be plotted.
%OUTPUTS
%   none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

default_options.draw_elements = 0;
default_options.element_edge_color = [1,1,1] * 0.95;
default_options.element_face_colour = [1,1,1] * 0.75;
default_options.mesh_edge_color = 'k';
default_options.draw_mesh_edges = 1;
default_options.node_sets_to_plot = [];
default_options.scale_factor = [];
default_options.el_mat_i = ones(size(els, 1), 1);
default_options.el_abs_i = zeros(size(els, 1), 1);
default_options.show_abs = 1;
default_options.matl_cols = [];
default_options.interface_el_col = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display_options = fn_set_default_fields(display_options, default_options);

if isempty(display_options.matl_cols)
    no_matls = numel(unique(display_options.el_mat_i));
    display_options.matl_cols = linspace(0.5, 0.75, no_matls)' * [1, 1, 1];
end

%Draw elements always to create the patch object. Then turn edges on or off
%etc to get desired effects
hold on;
els(els == 0) = NaN;
display_options.el_mat_i(display_options.el_mat_i == 0) = 1; %needed to prevent errors for elements with no material
base_cdata = permute(display_options.matl_cols(display_options.el_mat_i, :), [1, 3, 2]);
h = patch('Faces', els, 'Vertices', nds, 'CData', base_cdata, 'FaceColor', 'flat');

if display_options.interface_el_col
    [i, ~] = find(els == 0 | isnan(els));
    x = [nds(els(i, 1), 1), nds(els(i, 2), 1)]';
    y = [nds(els(i, 1), 2), nds(els(i, 2), 2)]';
    plot(x, y, display_options.interface_el_col);
end

if display_options.show_abs
    set(h, 'CData', base_cdata .* (1 - display_options.el_abs_i / 2))
end

if display_options.draw_elements
    set(h, 'EdgeColor', display_options.element_edge_color);
else
    set(h, 'EdgeColor', 'none');
end

if display_options.draw_mesh_edges
    %find edges that only occur once (i.e. they are the free edges)
    free_ed = fn_find_free_edges(els);
    %plot them
    hold on;
    plot(reshape(nds(free_ed, 1), size(free_ed))', reshape(nds(free_ed, 2), size(free_ed))', display_options.mesh_edge_color);
end

if ~isempty(display_options.node_sets_to_plot)
    hold on;
    for ii = 1:length(display_options.node_sets_to_plot)
        plot(nds(display_options.node_sets_to_plot(ii).nd, 1), nds(display_options.node_sets_to_plot(ii).nd, 2), display_options.node_sets_to_plot(ii).col);
    end
end

if display_options.draw_axes
    % Add axis lines
    ax = gca;
    box(ax, 'on');
    
    xlim([0, 1.5e-3]);

    aspect_ratio = [1, 4, 1]; % Adjust the second value to change the aspect ratio
    pbaspect(aspect_ratio);
  
    view(-90, 90); % Rotate around the x-axis by 90 degrees, then flip the view vertically

    set(gca,'FontName','Times','FontSize',8);

    xticks([0, 0.5e-3, 1e-3, 1.5e-3]); % Adjust tick positions as needed
    xticklabels({'0', '0.5', '1', '1.5'}); % Adjust tick labels as needed

    yticks([0, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3]); % Adjust tick positions as needed
    yticklabels({'6', '5', '4', '3', '2', '1', '0'}); % Adjust tick labels as needed

    xlabel('y-coordinate (mm)', 'Interpreter', 'latex');
    ylabel('x-coordinate (mm)', 'Interpreter', 'latex');
    
    % % Add text at specific coordinates
    % x_text = 0.75e-3; % X-coordinate
    % y_text = 2e-3; % Y-coordinate
    % text(x_text, y_text, 'G2','FontName','Times','FontSize',10, 'Interpreter', 'latex');
    % 
    % % Add text at specific coordinates
    % x_text = 0.75e-3; % X-coordinate
    % y_text = 5e-3; % Y-coordinate
    % text(x_text, y_text, 'G1','FontName','Times','FontSize',10, 'Interpreter', 'latex');

    hold on;
    yline(3.5e-3, ':','LineWidth', 1.5); % xline function automatically adjusts the y-limits
    hold off;

end
% 
% xlim([0, 3e-5]); % Set x-axis limits from 0 to 5e-5
% ylim([5e-4, 5.4e-4]);   % Set y-axis limits as needed

axis equal;
axis off;
end




function new_struct = fn_set_default_fields(old_struct, default_struct)
%USAGE
%	new_struct = fn_set_default_fields(old_struct, default_struct);
%SUMMARY
%	Use to add default fields and values to a structured variable, such as
%	options for a function.
%AUTHOR
%	Paul Wilcox (Dec 2003)
%INPUTS
%	old_struct - original structured variable
%	default_struct - structured variable containing default fields and
%	values
%OUTPUTS
%	new_struct - updated structured variable. All existing fields and values
%	in old_struct will be preserved, but any fields found in default_struct
%	and their values will be added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_struct = old_struct;
default_fieldnames = fieldnames(default_struct);
for ii=1:length(default_fieldnames)
	if ~isfield(new_struct, default_fieldnames{ii})
		new_struct = setfield(new_struct, default_fieldnames{ii}, getfield(default_struct, default_fieldnames{ii}));
    end
end
end
