import re

def write_graph_lines_and_get_ref(f_in, f_out):
    for line in f_in:
        if not line.startswith('W'):
            f_out.write(line)
        else:
            return line

pattern_walk_step = re.compile(r'([><])(\d+)')

def iter_walk(walk):
    for match in pattern_walk_step.finditer(walk):
        orientation = match.group(1)
        segment = int(match.group(2))
        yield orientation, segment



print('Enter file name:')
in_filename = input()
out_filename = "modified_" + in_filename

with open(in_filename, 'r') as f_in, open(out_filename, 'w') as f_out:

    ref = write_graph_lines_and_get_ref(f_in, f_out)
    f_out.write(ref)
    ref_steps = [(m[0], int(m[1:])) for m in re.findall(r'[<>]\d+', ref)]
    del ref

    for line in f_in:

        new_walk_steps = [] 
        buffer_walk_steps = []
        description, walk = line.rstrip('\r\n').rsplit("\t", 1)

        f_out.write(description + "\t")

        k = 0
        j = 0
        start_k = None
        start_j = None


        iterator_walk = iter_walk(walk)

        try:
            step_walk = next(iterator_walk)
        except StopIteration:
            f_out.write("\r\n")
            continue;
        
        while k < len(ref_steps):
            
            
            if ref_steps[k][1] == step_walk[1]:
                
                if start_k is not None and (k > start_k + 1 or j > start_j + 1):
                    new_walk_steps.extend(ref_steps[start_k + 1 : k])
                    buffer_walk_steps = []
                
                new_walk_steps.append(step_walk)

                start_k = k
                start_j = j
            
                k += 1
                try:
                    step_walk = next(iterator_walk)
                except StopIteration:
                    step_walk = None
                    break
                j += 1
                
            elif ref_steps[k][1] < step_walk[1]:
                k += 1
            else:
                if start_k is None:
                    new_walk_steps.append(step_walk)
                else:
                    buffer_walk_steps.append(step_walk)
                try:
                    step_walk = next(iterator_walk)
                except StopIteration:
                    step_walk = None
                    break
                j += 1

        if buffer_walk_steps:
            new_walk_steps.extend(buffer_walk_steps)

        if step_walk is not None:
            new_walk_steps.append(step_walk)

        new_walk = "".join(f"{orientation}{segment}" for orientation, segment in new_walk_steps) + "".join(f"{orientation}{segment}" for orientation, segment in iterator_walk)
        f_out.write(new_walk+"\r\n")