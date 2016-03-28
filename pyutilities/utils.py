from functools import wraps

def make_iterable(argn=0, exception=None):
    """Check to make sure an argument is an iterator.
    Default argument is args[0]
    Add exceptions that you wish to enclose in an iterator by providing an
    iterator of their type.
    Ironically, exceptions must be encapsulated in an iterator"""
    def decorate(func):
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            """Make a function iterable if it is not by putting it in a list.
            Note also converts a string to a list with the string inside"""
            maybe_iter = args[argn]
            if not hasattr(maybe_iter, '__iter__') or isinstance(maybe_iter, exception):
                maybe_iter = [maybe_iter]
            args = list(args)
            args[argn] = maybe_iter
            return func(*args, **kwargs)
        return func_wrapper
    return decorate
