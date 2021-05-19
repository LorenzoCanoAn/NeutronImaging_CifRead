function data = parameters_for_structure_factor(path_cif, path_vesta)
% Esta función acepta el path a un fichero de tipo .cif, y retorna una
% estructura con los siguientes componentes:
% - cell: contiente la información sobre la celdilla unidad
% - atom_site: contiente información en forma de tabla sobre el
% posicionamiento de los atomos en la red
% - atom_site_aniso: Es una tabla de 0x0 de tamaño si no se provee esta
% información en el .cif. Si la contiene, se expresa de la misma forma que
% atom_site.

%% Leer los archivos y meter todo el contenido en variables.
if isfile(path_cif)
    array_of_cif = file2array(path_cif);
    data.cif = process_cif_array(array_of_cif);
else
    data.cif = [];
end

if isfile(path_vesta)
    array_of_vesta = file2array(path_vesta);
    data.vesta = process_vesta_array(array_of_vesta);
else
    data.vesta = [];
end


end

%% Función que coje el archivo y lo almacena linea por linea en un array de celdas
function o_array = file2array(i_filePath)
% inputs:
%   i_filePath: string que contiene el path a el archivo por leer
% outputs:
%   o_array: array de deldas en el que cada celda es una linea del fichero
%            además, lo parentesis que indican la confianza en los
%            decimales se han eliminado para evitar problemas con la
%            conversión de string a numero posteriormente
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
%% Procesamiento de la información del archivo cif
% Recomendaría que si lo que se quiere es ampliar la función, añadir
% nuevas funciones aquí, en vez de modificar las ya existentes.
data.cell = get_cell_info(array_of_file);
data.atom_site = get_atom_type_info(array_of_file);
data.atom_site_aniso = get_atom_site_aniso_info(array_of_file);
data.sym_op = get_symetry_operations(array_of_file);
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
atom_site_aniso_info = cell(1,1);
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

function  vesta_info = process_vesta_array(array_of_file)
n_line = 0;
cell_of_data = cell(1,1);
for linea = array_of_file'
    n_line = n_line + 1;
    % Obtener los nombres de la primera linea para saber que significan
    % las lineas posteriores
    if n_line == 1
        nombres = split(linea)';
        void_els = nombres == "";
        nombres(void_els)=[];
        % Para las siguientes lineas, separar las lineas en función de los
        % espacios, y convertir cada elemento en un numero. Se da por hecho
        % que todos los element
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
function symOp = get_symetry_operations(array_of_file)
% inputs:
%   array_of_file: cell array that contains in each cell a line of the
%   original .cif file
% outputs:
%   - symOp: struct containing two fields:
%       - rawData: a cell array that contains on each line the string
%       representation of the symetri operations as written in the file
%       - matrices: a cell array that contains on each line the
%       transformation matrix corresponding to the symetry operation writen
%       in the same line in rawData.
n_line = 0;
inside_correct_loop = false;
symOp.rawData = cell(1,1);
symOp.matrices = cell(1,1);
for line = array_of_file' % loop over the file
    if contains(line, "_space_group_symop") % check when we enter the symetry operations part of the file
        inside_correct_loop = true;
    end
    if inside_correct_loop % if we are inside the symetry operations part
        % check if we have finished, if so, exit the function
        if contains(line,"loop_")
            break
        end
        % check if we are finally on a symetry operation. If that is
        % the case, the splitting of the line should give an array of
        % more than one element
        sp_line = split(line);
        if max(size(sp_line)) > 3
            n_line = n_line + 1;
            symOp.rawData{n_line,1} = line;
            [T_matrix, id] = symOp_name_to_matrix(line);
            symOp.matrices{n_line,1} = T_matrix;
            symOp.matrices{n_line,2} = id;
        end
        
        
    end
end
end

function [T_matrix, id]= symOp_name_to_matrix(string_rep)
string_rep = strrep(string_rep,"'","");
string_rep = strrep(string_rep,",","");
s_str = split(string_rep);
id = str2num(s_str(1));
s_str(1) = [];

% at this point, only 3 elements should be in the splitted string, each
% regarding, regarding x,y and z.
if max(size(s_str)) ~= 3
    error(strcat("There is a problem with the symop with id: ",num2srt(id)));
end

T_matrix = zeros(3,4);
list_coord = ["x";"y";"z"];
%% Rotation matrix
for i = 1:3 % loop over the three elements of s_str
    info = s_str(i);
    info_c = convertStringsToChars(info);
    for j = 1:3 % loop over list_coord
        k = strfind(info,list_coord(j));
        if max(size(k)) > 0 % if the coordinate is present in this string
            if k > 1
                if info_c(k-1)=="-" % if the coordinate is negative
                    T_matrix(i,j) = -1;
                else
                    T_matrix(i,j) = 1;
                end
            else
                T_matrix(i,j) = 1;
            end
        end
    end
end
%% Translation vector
list_numbers = "1234567890/";
for i = 1:3 % loop over the three elements of s_str
    only_num = [];
    count = 0;
    begin_index = -1;
    char_array = convertStringsToChars(s_str(i));
    for c = char_array
        count = count + 1;
        if contains(list_numbers,c)
            if begin_index < 0
                begin_index = count;
            end
            only_num = cat(2,only_num,c);
        end
    end
    
    
    
    if max(size(only_num)) == 0
        continue
    else
        sign = 1;
        if max(size(char_array)) > 1
            if char_array(begin_index-1) == "-"
                sign = -1;
            end
        end
        if contains(only_num,"/")
            splitted = split(only_num,"/");
            if max(size(splitted)) == 2
                T_matrix(i,4) = str2double(splitted{1})/str2double(splitted{2})*sign;
            else
                error(strcat("There is a problem with the symop with id: ",num2srt(id)));
            end
        else
            T_matrix(i,4) = str2double(only_num)*sign;
        end
    end
end

end