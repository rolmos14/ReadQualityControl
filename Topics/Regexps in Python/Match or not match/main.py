import re


def matched(template, string):
    return re.match(template, string) is not None
