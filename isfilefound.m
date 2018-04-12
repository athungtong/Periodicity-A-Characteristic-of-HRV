function flag=isfilefound(filepath)
flag=1;
fid=fopen(filepath); 
if fid==-1
    display(['filepath: ' filepath ' not found.'])
    flag=0;
else
    fclose(fid);
end
