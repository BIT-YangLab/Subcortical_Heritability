function results = readlines(filename)
results = [];
fid = fopen(filename);

while ~feof(fid)
    str = fgetl(fid);
    results = [results; str];
end

fclose(fid);