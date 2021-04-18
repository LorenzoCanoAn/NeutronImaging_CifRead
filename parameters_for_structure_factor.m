function data = parameters_for_structure_factor(path_cif, path_vesta)
% Esta función acepta el path a un fichero de tipo .cif, y retorna una
% estructura con los siguientes componentes:
% - cell: contiente la información sobre la celdilla unidad
% - atom_site: contiente información en forma de tabla sobre el
% posicionamiento de los atomos en la red
% - atom_site_aniso: Es una tabla de 0x0 de tamaño si no se provee esta
% información en el .cif. Si la contiene, se expresa de la misma forma que
% atom_site.

%% Get content of files
array_of_cif = file2array(path_cif);
array_of_vesta = file2array(path_vesta);
data.cif = process_cif_array(array_of_cif);
data.vesta = process_vesta_array(array_of_vesta);
end

%% Function to put the file into an string array
function o_array = file2array(i_filePath)
% The input is the path to the file to read. The output is an array where
% each row of the array is a string that corresponds to a line of the file.
fileID = fopen(i_filePath);
o_array = [];
while true
    line = fgetl(fileID);
    if line == -1 % If the end of the file has been reached
        break
    end
    s_line = string(line);
    s_line = strrep(s_line,"(","");
    s_line = strrep(s_line,")","");
    o_array = cat(1,o_array,s_line);
end
fclose(fileID);
end

%% Function to process the data array of the file
function data = process_cif_array(array_of_file)
    data.cell = get_cell_info(array_of_file);
    data.atom_site = get_atom_type_info(array_of_file);
    data.atom_site_aniso = get_atom_site_aniso_info(array_of_file);
end

function cell_info = get_cell_info(array_of_file)
for line = array_of_file'
    splitted_line = split(line);
    switch splitted_line(1)
        case "_cell_length_a"
            cell_info.a = str2double(splitted_line(2));
        case "_cell_length_b"
            cell_info.b = str2double(splitted_line(2));
        case "_cell_length_c"
            cell_info.c = str2double(splitted_line(2));
        case "_cell_angle_alpha"
            cell_info.alpha = str2double(splitted_line(2));
        case "_cell_angle_beta"
            cell_info.beta = str2double(splitted_line(2));
        case "_cell_angle_gamma"
            cell_info.gamma = str2double(splitted_line(2));
        case "_cell_volume"
            cell_info.vol = str2double(splitted_line(2));
        case "_cell_formula_units_Z"
            cell_info.z_units = str2double(splitted_line(2));
        case "_space_group_name_H-M_alt"
            cell_info.sg_name = join(splitted_line(2:end)," ");
        case "_space_group_IT_number"
            cell_info.IT_num= str2double(splitted_line(2));
    end
end
end

function atom_type_info = get_atom_type_info(array_of_file)
    inside_loop = false;
    atom_type_info = cell(1,1);
    n_fields = 0;
    n_elements = 0;
    for line = array_of_file'
        if contains(line,"_atom_site_") && ~contains(line,"_atom_site_aniso_")
            n_fields = n_fields + 1;
            inside_loop = true;
            name = strrep(line,"_atom_site_","");
            atom_type_info{1,n_fields} = name;
        else
            if inside_loop
                char_line = char(line);
                if char_line(1)=="_" || contains(line,"loop_")
                    % if the loop over the atom site info has 
                    % ended, break
                    break
                elseif char_line(1)=="#"
                    continue
                else
                    n_elements = n_elements + 1;
                    n_field = 0;
                    split_line = split(line);
                    for element = split_line'
                        n_field = n_field + 1;
                        number_element = str2double(element);
                        if isnan(number_element)
                            atom_type_info{n_elements+1,n_field}=element;
                        else
                            atom_type_info{n_elements+1,n_field}=number_element;
                        end
                    end
                end
            end
        end
    end
    c = atom_type_info(2:end,:);
    names = [atom_type_info{1,:}];
    atom_type_info = cell2table(c,'VariableNames',names);
end

function atom_site_aniso_info = get_atom_site_aniso_info(array_of_file)
    inside_loop = false;
    atom_site_aniso_info = {};
    n_fields = 0;
    n_elements = 0;
    for line = array_of_file'
        if contains(line,"_atom_site_aniso")
            n_fields = n_fields + 1;
            inside_loop = true;
            name = strrep(line,"_atom_site_aniso_","");
            atom_site_aniso_info{1,n_fields} = name;
        else
            if inside_loop
                char_line = char(line);
                if char_line(1)=="_" || contains(line,"loop_")
                    % if the loop over the atom site info has 
                    % ended, break
                    break
                elseif char_line(1)=="#"
                    continue
                else
                    n_elements = n_elements + 1;
                    n_field = 0;
                    split_line = split(line);
                    for element = split_line'
                        n_field = n_field + 1;
                        number_element = str2double(element);
                        if isnan(number_element)
                            atom_site_aniso_info{n_elements+1,n_field}=element;
                        else
                            atom_site_aniso_info{n_elements+1,n_field}=number_element;
                        end
                    end
                end
            end
        end
    end
    if isempty(atom_site_aniso_info)
        atom_site_aniso_info = table;
    else
        c = atom_site_aniso_info(2:end,:);
        names = [atom_site_aniso_info{1,:}];
        atom_site_aniso_info = cell2table(c,'VariableNames',names);
    end
end

function  vesta_info = process_vesta_array(data_array)
    n_line = 0;
    cell_of_data = cell(1,1);
    for linea = data_array'
        n_line = n_line + 1;
        if n_line == 1
            nombres = split(linea)';
            void_els = nombres == "";
            nombres(void_els)=[];
        else
            elementos = split(linea)';
            void_els = elementos == "";
            elementos(void_els)=[];
            n_el = 0;
            for elemento = elementos
                element_num = str2double(elemento);
                n_el = n_el+1;
                if isnan(element_num)
                    cell_of_data{n_line-1,n_el} = elemento;
                else
                    cell_of_data{n_line-1,n_el} = element_num; 
                end
            end
        end
    end
    vesta_info = cell2table(cell_of_data,'VariableNames',nombres);
end