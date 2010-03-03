% test for .ov file read/write

name = 'OVTable_2Handle';

[vertex,face] = read_ov([name '.ov']);

dummy = write_ov([name '-tmp.ov'], vertex,face);

[vertex1,face1] = read_ov([name '-tmp.ov']);