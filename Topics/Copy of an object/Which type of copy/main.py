import copy


def detect_copy():
    obj = [[1]]
    obj_copy = copying_machine(obj)
    if obj[0] is obj_copy[0] and obj[0][0] is obj_copy[0][0]:
        return "shallow copy"
    return "deep copy"
