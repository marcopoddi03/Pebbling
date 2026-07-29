import re
import sys

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


if len(sys.argv) > 1:
    in_filename = sys.argv[1]
else:
    print('Enter file name:')
    in_filename = input()

out_filename = "revised_" + in_filename

with open(in_filename, 'r') as f_in, open(out_filename, 'w') as f_out:

    ref = write_graph_lines_and_get_ref(f_in, f_out)
    f_out.write(ref)
    ref_walk = ref.rstrip('\r\n').rsplit("\t", 1)[1]

    for line in f_in:

        buffer_walk_steps = []
        buffer_ref_steps = []
        description, walk = line.rstrip('\r\n').rsplit("\t", 1)

        f_out.write(description + "\t")

        k = 0
        j = 0
        start_k = None
        start_j = None


        iterator_walk = iter_walk(walk)
        iterator_ref = iter_walk(ref_walk)

        try:
            step_walk = next(iterator_walk)
        except StopIteration:
            f_out.write("\r\n")
            continue;

        try:
            step_ref = next(iterator_ref)
        except StopIteration:
            break
        
        while True:
            
            
            if step_ref[1] == step_walk[1]:
                
                if start_k is not None and (k > start_k + 1 or j > start_j + 1):
                    for m in buffer_ref_steps:
                        f_out.write(f"{m[0]}{m[1]}")
                    buffer_walk_steps = []
                
                f_out.write(f"{step_walk[0]}{step_walk[1]}")

                buffer_ref_steps = []

                start_k = k
                start_j = j
            
                k += 1
                try:
                    step_walk = next(iterator_walk)
                except StopIteration:
                    step_walk = None
                    break
                j += 1
                try:
                    step_ref = next(iterator_ref)
                except StopIteration:
                    break
                
            elif step_ref[1] < step_walk[1]:
                if start_k is not None:
                    buffer_ref_steps.append(step_ref)
                try:
                    step_ref = next(iterator_ref)
                except StopIteration:
                    break
                k += 1
            else:
                if start_k is None:
                    f_out.write(f"{step_walk[0]}{step_walk[1]}")

                else:
                    buffer_walk_steps.append(step_walk)
                try:
                    step_walk = next(iterator_walk)
                except StopIteration:
                    step_walk = None
                    break
                j += 1

        if buffer_walk_steps:
            for m in buffer_walk_steps:
                f_out.write(f"{m[0]}{m[1]}")
        if step_walk is not None:
            f_out.write(f"{step_walk[0]}{step_walk[1]}")

        new_walk = "".join(f"{orientation}{segment}" for orientation, segment in iterator_walk)
        f_out.write(new_walk+"\r\n")