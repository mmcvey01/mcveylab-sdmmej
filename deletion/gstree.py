# -*- coding: utf-8 -*-
"""
Generalized suffix trees

These are compact suffix trees, i.e., (non-root) internal nodes always
branch.  Construction uses Ukkonen's on-line linear-time algorithm.
"""

__credits__ = """\
References
==========

Baker, Brenda S. "Parameterized Duplication in Strings: Algorithms and
an Application to Software Maintenance", _SIAM J. Comput._ 26 (1997):
1343-1362.

Baker, Brenda S. "On finding duplication and near-duplication in large
software systems." _Proceedings of the 2nd Working Conference on
Reverse Engineering, 1995_: 86-95.

Giegerich, Robert and Stefan Kurtz. "From Ukkonen to McCreight and
Weiner: A Unifying View of Linear-Time Suffix Tree Construction."
_Algorithmica_ 19 (1997): 331-353.

Gusfield, Dan. _Algorithms on Strings, Trees, and Sequences_.
Cambridge: Cambridge University Press, 1997.

Maaß, Moritz G. "Linear Bidirectional On-Line Construction of Affix
Trees." _Proc. Combinatorial Pattern Matching 2000, Lecture Notes in
Computer Science_ 1848 (2000): 320-334.

Stoye, Jens. "Affix Trees." Technical Report 2000-04. Berlin:
Universität Bielefeld, Technische Fakultät, 2000.

Ukkonen, Esko. "On-line construction of suffix trees." _Algorithmica_
14 (1995): 249-260.
"""

__copyright__ = """Copyright (C) 2009-2011 Tom Yu <tlyu@mit.edu>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holders nor the names of their
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

__author__ = 'Tom Yu <tlyu@mit.edu>'

__all__ = 'GSTree'

_debug = False

# Use booleans, so that the "not" operator works to flip views.
PREFIX = False
SUFFIX = True

from UserDict import DictMixin
from itertools import *

# "LList" for linked lists; "list" for Python builtin lists.  This
# improves maxpairs performance for strings of more than a few
# thousand characters.

# from llist import LList
_list = list

# Utility functions

def _classstr(obj):
    cls = obj.__class__
    return "%s.%s" % (cls.__module__, cls.__name__)

def _pathlen(parents):
    length = 0
    for edge in parents:
        length += len(edge)
    return length

def _slice2ndx(sl):
    L = [str(x) if x is not None else '' for x in sl.start, sl.stop]
    if sl.step == -1:
        L.append(str(sl.step))
    return ':'.join(L)

def _view2step(view, __step=(-1, 1)): return __step[view]
def _step2view(step, __view=(PREFIX, None, SUFFIX)): return __view[step+1]

class TreeStr(object):
    """A string (or other string-like sequence object) in a tree.

    It knows what its tree, key, and terminator are.
    """
    __slots__ = ('tree', 'key', 's', 'term', 'init', 'sbound', 'pbound')

    def __init__(self, tree, key, s):
        self.tree = tree
        self.key = key
        self.s = s
        self.term = tree._key2term(key)
        self.init = tree._key2init(key)
        self.sbound = 0         # asp.stop
        self.pbound = -len(s) - 1 # app.stop

    def __repr__(self):
        return ("<%s %s tree=%s key=%s s=%s>" %
                (_classstr(self), hex(id(self)), hex(id(self.tree)),
                 self.key, self.s))

    def __setitem__(self, k, v):
        raise TypeError("item setting not supported (yet)")
    def __iadd__(self, v):
        raise TypeError("in-place append not supported (yet)")

    def __getattr__(self, attr):
        return getattr(self.s, attr)

    def __str__(self):          return str(self.s)
    def __len__(self):          return len(self.s)

    # __getitem__ doesn't properly raise IndexError on out-of-bounds
    # indices because of boundary pseudo-character behavior.  Define
    # an __iter__ method in case some caller wants to iterate over a
    # TreeStr, to avoid infinite loops.
    def __iter__(self):         return iter(self.s)

    def __getitem__(self, k):
        try:
            return self.s[k]
        except IndexError:
            slen = len(self.s)
            if k >= slen:
                # positive indices past the end get "$" terminator char
                return self.term
            elif k < -slen:
                # negative indices past the beginning get "^" initiator char
                return self.init
            else:
                raise


class Word(object):
    """A word in a tree.  Can correspond to an edge or a node reference.

    The word is specified by a key (indicating which string) and slice
    indices into that string.  A slice stop of None indicates that the
    word extends to the end of the string (an "open edge").
    """
    # Unfortunately, can't subclass "slice" (yet?), but that's
    # probably not a bad thing beause slice objects are immutable.
    __slots__ = ('tstr', 'start', 'stop', 'step')

    def __init__(self, tstr, start, stop, step):
        # if not isinstance(tstr, TreeStr):
        #     raise TypeError("'%s' object is not a TreeStr" % type(tstr))
        self.tstr = tstr
        self.start = start
        self.stop = stop
        self.step = step

    def __len__(self):
        slen = len(self.tstr.s)
        start, stop, step = self.start, self.stop, self.step
        if start is None:
            start = 0
        elif start < 0:
            start += slen + 1
        if stop is None:
            stop = slen
        elif stop < 0:
            stop += slen + 1
        stop = min(max(stop, -1), slen)
        start = min(max(start, -1), slen)
        return max((stop - start) * step, 0)

    def __str__(self):
        return str(self.word())

    def __repr__(self):
        return "<%s word=%s>" % (_classstr(self), self.strword())

    @property
    def adjstart(self):
        start = self.start
        if start is not None:
            return start
        step = self.step
        if step == -1:
            start = self.tstr.sbound
        else:
            start = self.tstr.pbound
        return start + (len(self.tstr.s) + 1) * step

    @property
    def adjstop(self):
        stop = self.stop
        if stop is not None:
            return stop
        step = self.step
        if step == -1:
            stop = self.tstr.pbound
        else:
            stop = self.tstr.sbound
        return stop

    @property
    def adjslice(self):
        start, stop, step = self.adjstart, self.adjstop, self.step
        if start < 0 and step == 1:
            start = 0
        elif start >=0 and step == -1:
            start = -1
        return slice(start, stop, step)

    @property
    def adjlen(self):
        return (self.adjstop - self.adjstart) * self.step

    def isopen(self):
        """Whether this word corresponds to an open edge."""
        return self.stop is None

    @property
    def letter(self):
        tstr, start, step = self.tstr, self.start, self.step
        if start is not None:
            return tstr[start]

        slen = len(tstr.s)
        if step == -1:
            start = self.tstr.sbound
            if start <= slen:
                start -= slen + 1
        else:
            start = self.tstr.pbound
            if start >= -slen - 1:
                start += slen + 1

        return tstr[start]

    def word(self):
        """The actual word, as a sequence.

        This is not necessarily a string if the tree is built from a
        non-string sequence type."""
        return self.tstr[self.adjslice]

    def strword(self):
        """Key, indices, and substring of this word.

        This is primarily a debugging aid.  Takes the form

            'key[start:stop:step]substring'

        where "start", "stop", and "step" are as in slice notation.
        """
        sl = slice(self.start, self.stop, self.step)
        return ("%s[%s]%s" %
                (self.tstr.key, _slice2ndx(sl), str(self.word())))

    def rev(self):
        tstr = self.tstr
        slen = len(tstr)
        start, stop, step = self.start, self.stop, self.step
        step *= -1
        if start is not None:   start += (slen + 1) * step
        if stop is not None:    stop += (slen + 1) * step
        return Word(tstr, stop, start, step)


class _Cwrap(object):
    def __init__(self, d):
        self.__d = d
    def __getattr__(self, attr):
        return self.__d.get(attr)


class Node(object):
    """An explicit node in the suffix tree.

    May be viewed as either a suffix node or a prefix node.

    A Node has a label that is a Word corresponding to the path
    leading to that Node from the root.  A leaf node has an open Word
    as its label.  The root node has no label.  The step value of the
    Word is +1 for the suffix view and -1 for the prefix view.

    The Word labeling a Node is reversed in its dual view.  The key
    and indices of the Word are probably not unique for an inner node,
    but may depend on the order in which strings are inserted into the
    tree.

    Edges are implied by the node labels of the nodes they join.

    Equivalences:

        prefix          suffix
        ---------------------------
        child           rchild
        rchild          child
        prefix link     parent
        parent          suffix link

    The children of a node are always have the same view as their
    parents.  The suffix (or prefix) link target of a node has the
    same view as the node.  This means that the parent of a node is
    has the opposite view compared to that node.

    Special case representations for node labels:

        * Full string leaf: [::1]    [::-1]
        * Terminator leaf:  [len::1] [:-1:-1]
        * Initiator leaf:   [:0:1]   [-len-1::-1]
    """

    __slots__ = ('rev', 'child', 'link', 'loc', 'zpath',
                 'lchr', 'lchrlists')

    def __init__(self, _node=None, _loc=None, *args, **kwargs):
        self.child = {}
        self.link = None
        self.zpath = None
        if _node is not None:
            self.rev = _node
            try:
                self.loc = _node.loc.rev()
            except AttributeError:
                self.loc = None
        else:
            self.loc = _loc
            self.rev = self._newrev()

    def _newrev(self):
        raise TypeError('Node cannot be instantiated directly; '
                        'use SNode or PNode.')

    def __repr__(self):
        loc, parent, link = self.loc, self.parent, self.link
        if loc is not None:     loc = loc.strword()
        if parent is not None:  parent = hex(id(parent))
        if link is not None:    link = hex(id(link))
        return ("<%s %s loc=%s child=%s rchild=%s parent=%s link=%s>" %
                (_classstr(self), hex(id(self)), loc,
                 self.child.keys(), self.rchild.keys(),
                 parent, link))

    def __invert__(self):
        """~node -> parent of node -- for internal/debugging use only
        """
        return self.parent.rev

    @property
    def parent(self):           return self.rev.link
    @property
    def rchild(self):           return self.rev.child

    @parent.setter
    def parent(self, val):      self.rev.link = val

    @property
    def _sortkey(self):
        try:
            return self.loc.tstr.tree._sortkey
        except AttributeError:
            # Prefix view of root has no loc, so get _sortkey from the
            # suffix view.
            return self.rev._sortkey

    @property
    def _(self):
        """Hack to make it easier to interactively walk the tree.

        t._.a._.b._.c <-> t.child['a'].child['b'].child['c']
        """
        return _Cwrap(self.child)

    def _add_leaf(self, tstr, i):
        """Add a new leaf node as a child.

        The path to the new leaf is an open edge.
        """
        step = _view2step(self._view)
        nleaf = self._newnode(Word(tstr, i, None, step))
        self._add_child(nleaf, collapse=False)
        return nleaf

    def _add_child(self, child, collapse=True):
        if collapse and len(child.child) == 1:
            self._add_zptail(child)
        # elif collapse and self.parent and self.parent.rev.zpath:
        #     self._add_zphead(child)
        else:
            self._add_disowned(child)
            self._adopt(child)

    def _adopt(self, child):
        child.parent = self.rev

    def _add_disowned(self, child):
        """Add a child, but don't update the child's parent pointer.
        """
        echr = self._edgechr(child)
        if echr in self.child:
            raise ValueError("%s already has a child for edgechr '%s'" %
                             (self, echr))
        self.child[echr] = child

    def _add_zphead(self, child):
        """Extend the head of the ZPath."""
        self._add_disowned(child)
        self._adopt(child)
        zpath = self.parent.rev.path
        zpath.addhead(child)

    def _add_zptail(self, child):
        """Extend the tail of the ZPath."""
        # Typically, this happens when creating a suffix link from one
        # node (A) to another (B), and takes the form of
        # B.rev._add_child(A.rev) -- really adding a child link from B
        # to A in the dual tree.
        self._adopt(child)
        grandchild = child.onlychild()
        zpath = grandchild.parent.rev.zpath
        if not zpath:
            child.zpath = zpath = ZPath(child)
        else:
            child.zpath = zpath
            zpath.addtail(child)

        zphead = zpath.head
        self._add_disowned(zphead)

    def _delchild(self, child):
        """Delete a child, but don't change its parent pointer."""
        echr = self._edgechr(child)
        del self.child[echr]

    def _edgechr(self, child):
        """First character along the edge from self to child.
        """
        s = child.loc.tstr
        start = child.loc.adjstart
        myloc = self.loc
        if myloc is None:
            return child.loc.letter
        else:
            i = start + myloc.adjlen * myloc.step
            return s[i]

    def _edgelen(self, child):
        """Length of the edge from self to child.
        """
        myloc = self.loc
        childloc = child.loc
        if myloc is None:
            return childloc.adjlen
        else:
            return childloc.adjlen - myloc.adjlen

    def _edgeword(self, child):
        """Word corresponding to the edge from self to child.
        """
        myloc = self.loc
        childloc = child.loc
        if myloc is None:
            return childloc
        else:
            step = childloc.step
            start = childloc.adjstart + myloc.adjlen * step
            stop = childloc.stop
            return Word(childloc.tstr, start, stop, step)

    def _split(self, child, wlen, collapse=False, recurse=False):
        """Split edge between self and child at wlen."""
        if _debug:
            print "_split:\n%s\n%s @ %s" % (self, child, wlen)
        if not recurse and child.parent.rev.zpath:
            return self._zsplit(child, wlen, collapse)

        nn = self._split_nn(child, wlen)

        self._delchild(child)

        nn._add_child(child, collapse=False)
        self._add_child(nn, collapse)
        return nn

    def _split_nn(self, child, wlen):
        """Create new node for insertion into split point.

        The edge to be split is the edge self->child.
        """
        step = _view2step(self._view)
        childloc = child.loc
        loc = self.loc
        stop = childloc.adjstart
        if loc is None:
            stop += wlen * step
        else:
            stop += (wlen + loc.adjlen) * step
        start = childloc.start
        return self._newnode(Word(childloc.tstr, start, stop, step))

    def _split_add_leaf(self, child, wlen, tstr, i):
        nn = self._split(child, wlen)
        nleaf = nn._add_leaf(tstr, i)
        return nn, nleaf

    def _zsplit(self, child, wlen, collapse):
        """Split a collapsed edge at wlen."""
        zpath = child.parent.rev.zpath
        ztail, zhead = zpath[0], zpath[-1]
        tlen = self._edgelen(ztail)
        hlen = zhead._edgelen(child)
        if wlen < tlen:
            # Split tail edge.
            nn = self._split_nn(ztail, wlen)
            nn._adopt(ztail)
            if collapse:
                nn._add_disowned(ztail)
                self._adopt(nn)
                zpath.addtail(nn)
            else:
                nn._add_disowned(child)
                self._delchild(child)
                self._add_child(nn, collapse=False)
            return nn

        elif wlen > tlen + len(zpath) - 1:
            # Split head edge.
            wlen -= (tlen + len(zpath) - 1)
            nn = zhead._split_nn(child, wlen)
            nn._add_child(child, collapse=False)
            zhead._delchild(child)
            zhead._add_child(nn, collapse=False)
            if collapse:
                zpath.addhead(nn)
            else:
                self._delchild(child)
                self._add_disowned(nn)
            return nn

        else:
            i = wlen - tlen
            nn = zpath.split(i, collapse)
            if collapse:
                return nn

            self._delchild(child)
            self._add_disowned(nn)
            nn._delchild(nn.onlychild())
            nn._add_disowned(child)
            return nn

    def assertonechild(self):
        if len(self.child) != 1:
            raise ValueError("%s violates one-child policy" % self)

    def findchild(self, c):
        return self.child.get(c)

    def haschild(self, c):
        return c in self.child

    def isleaf(self):
        return not self.child

    def isroot(self):
        return self.loc is None

    def iterchildren(self):
        return self.child.itervalues()

    def onlychild(self):
        self.assertonechild()
        return next(self.child.itervalues())

    def sortedchildren(self):
        keys = sorted(self.child.keys(), key=self._sortkey)
        return (self.child[k] for k in keys)

    def printleaves(self, **kwargs):
        """Print the leaf nodes of this tree (or subtree).

        >>> t = GSTree('mississippi')
        >>> t.printleaves(hexids=False)
        0[11:] (len=0)
        0[1:2]i 0[11:] (len=1)
        0[1:2]i 0[8:]ppi (len=4)
        0[1:2]i 0[2:5]ssi 0[8:]ppi (len=7)
        0[1:2]i 0[2:5]ssi 0[5:]ssippi (len=10)
        0[:]mississippi (len=11)
        0[8:9]p 0[10:]i (len=2)
        0[8:9]p 0[9:]pi (len=3)
        0[2:3]s 0[4:5]i 0[8:]ppi (len=5)
        0[2:3]s 0[4:5]i 0[5:]ssippi (len=8)
        0[2:3]s 0[3:5]si 0[8:]ppi (len=6)
        0[2:3]s 0[3:5]si 0[5:]ssippi (len=9)
        """
        for x in IterLeafStrs(self, **kwargs):
            print x

    def printnodes(self, **kwargs):
        """Print all nodes of this tree (or subtree)."""
        for x in IterNodeStrs(self, **kwargs):
            print x

    def iterwhere(self, word):
        """Iterator over all locations of a word in strings of the tree.

        Locations are (key, index) tuples.
        """
        return IterWhere(self, word)

    def where(self, word):
        """List of all locations of a word in strings of the tree.

        Locations are (key, index) tuples.
        """
        return list(self.iterwhere(word))

    def itermaxpairs(self, *args, **kwargs):
        """Iterator over maximal pairs of this tree (or subtree).

        Keyword argument "minlen" (defaults to 0) is the minimum
        length of each repeated word.

        Keyword argument "allpairs" (defaults to False), if True,
        returns all pairs, not just maximal pairs.

        A maximal pair is a pair of occurrences of a substring such
        that extending each occurrence to the left or to the right
        will cause the pair to cease matching.

        Each maximal pair is a tuple

            ((key1, index1), (key2, index2), wordlen)

        where key1 and key2 are keys denoting strings in the tree,
        index1 and index2 are the indices (zero-based) into those
        strings where the word starts, and wordlen is the length of
        the word.

        >>> t = GSTree('mississippi')
        >>> sorted(t.itermaxpairs()) # doctest: +NORMALIZE_WHITESPACE
        [((0, 1), (0, 4), 4), ((0, 1), (0, 7), 1), ((0, 1), (0, 10), 1),
        ((0, 2), (0, 3), 1), ((0, 2), (0, 6), 1), ((0, 3), (0, 5), 1),
        ((0, 4), (0, 10), 1), ((0, 5), (0, 6), 1), ((0, 7), (0, 10), 1),
        ((0, 8), (0, 9), 1)]
        """
        return IterMaxPairs(self, *args, **kwargs)

    def maxpairs(self, *args, **kwargs):
        """List of maximal pairs of this tree (or subtree).

        See itermaxpairs() for details.
        """
        return list(self.itermaxpairs(*args, **kwargs))


class PNode(Node):
    __doc__ = Node.__doc__
    __slots__ = ()

    @staticmethod
    def _newnode(loc):          return PNode(_loc=loc)

    def _newrev(self):          return SNode(_node=self)

    _view = PREFIX

class SNode(Node):
    __doc__ = Node.__doc__
    __slots__ = ()

    @staticmethod
    def _newnode(loc):          return SNode(_loc=loc)

    def _newrev(self):          return PNode(_node=self)

    _view = SUFFIX


class ZPath(object):
    """A collapsed path.

    A collapsed path is a cluster of nodes that forms a linear
    subgraph and behaves as a single edge.
    """
    __slots__ = ('_nodes', 'head', 'tail')

    def __init__(self, node=None):
        self._nodes = []
        if node is None:
            return
        self.addtail(node)

    def __getitem__(self, i):   return self._nodes[i]
    def __len__(self):          return len(self._nodes)

    @property
    def head(self):     return self._nodes[-1].onlychild()
    @property
    def tail(self):     return self._nodes[0].parent

    def addtail(self, node):
        node.assertonechild()
        if len(self._nodes) > 1:
            self._nodes[0].zpath = None
        self._nodes.insert(0, node)
        node.zpath = self

    def addhead(self, node):
        node.assertonechild()
        if len(self._nodes) > 1:
            self._nodes[-1].zpath = None
        self._nodes.insert(len(self._nodes), node)
        node.zpath = self

    def split(self, i, collapse):
        n = self._nodes[i]
        if collapse:
            return n

        if len(self._nodes) == 1:
            n.zpath = None
            return n

        if i == len(self._nodes) - 1:
            del self._nodes[-1]
            self._nodes[-1].zpath = self
        elif i == 0:
            del self._nodes[0]
            self._nodes[0].zpath = self
        else:
            nzp = ZPath()
            nzp._nodes = self._nodes[i+1:]
            del self._nodes[i:]
            nzp._nodes[-1].zpath = nzp
            self._nodes[-1].zpath = self
            nzp._nodes[0].zpath = nzp

        n.zpath = None
        return n

class NodeRef(Word):
    """A reference pair for a (possibly implicit) node.

    A reference pair consists of a Node (an explicit node) and a
    (possibly empty) word.  A word is represented by a slice into a
    string, specifying a substring.

    Implicit nodes are not Node objects; they are "hidden" inside an
    Edge, and thus require an additional word for a full description.

    If the slice is "empty", the word is the empty string and this
    NodeRef is an explicit reference, i.e. to an actual Node object.
    """
    __slots__ = ('node', 'pathlen')

    def __init__(self, node, tstr, start, stop, step, pathlen=0):
        self.node = node
        self.pathlen = pathlen
        super(NodeRef, self).__init__(tstr, start, stop, step)

    def __repr__(self):
        return ("<%s node=%s word=%s>" %
                (_classstr(self), hex(id(self.node)), self.strword()))

    def isexplicit(self):
        """Whether this reference pair specifies an explicit node."""
        return self.start == self.stop

    def canonize(self):
        """Make this NodeRef canonical.

        If the path from the explicit node of this NodeRef to the
        (possibly implicit) node specified by this NodeRef traverses
        any explicit nodes, update self.node and self.start
        until that is no longer true, while continuing to refer to the
        same node.
        """
        s, node, start, step = self.tstr, self.node, self.start, self.step
        wlen = len(self)
        while True:
            if wlen <= 0: break
            c = s[start]
            child = node.findchild(c)
            if child.isleaf():
                break
            elen = node._edgelen(child)
            if elen > wlen:
                break
            start += elen * step
            wlen -= elen
            node = child

        self.node, self.start = node, start

    def followlink(self):
        """Follow the suffix or prefix link of this NodeRef.

        This may involve simulating the suffix link of an internal
        node.
        """
        implicit = not self.isexplicit()
        if self.node.link is not None:
            self.node = self.node.link
            self.pathlen -= 1
        else:
            # Root node doesn't have a link.
            if implicit:
                self.start += self.step
                self.pathlen -= 1
        if implicit:
            self.canonize()


class ActiveRef(NodeRef):
    """The reference pair of the active suffix (or prefix).

    This is the the reference pair of the longest nested suffix (or
    prefix) of the portion of the string that is already inserted.  A
    suffix (prefix) is nested if it occurs more than once.

    This contains the methods necessary to extend the tree one
    character at a time.

    "The node (u,v) representing the active suffix of t in cst(t) is
    the neuralgic point of the suffix tree." (Giegerich and Kurtz 1997)
    """

    def __init__(self, node, tstr, start, stop, step):
        super(ActiveRef, self).__init__(node, tstr, start, stop, step)
        self.activeleaf = None
        self.fullstr = None

    def update(self, i):
        """Update the active suffix (or prefix), given a position in a
        string in this tree.

        For now, this must be done in left-to-right order, and the
        reference pair must already be canonical.

        self.tstr[self.stop] is the character being inserted.

        self.tstr[self.stop-self.step] is the previous character that
        was inserted.

        self.pathlen is the total length of the active ref, which may
        be greater than len(self), which is the length of the word
        leading from the explicit node.
        """
        if _debug:
            print "\nupdate(%d): chr=%s" % (i, self.tstr[i])
        self._unlink_activeleaf()
        if self.step == -1:
            self.tstr.pbound -= 1
        else:
            self.tstr.sbound += 1
        # tracks the most recently added node, to ease suffix linking
        self.prevnode = None
        end = False
        self._split_fullstr()
        while not end:
            assert self.stop == i
            if self.isexplicit():
                end = self._update_explicit()
            else:
                end = self._update_implicit()
        self._term_activeleaf()

    def _link_prevnode(self, node):
        """Update the suffix link of the previously added node, i.e.,
        set the suffix link of prevnode to point to node."""

        prevnode = self.prevnode
        if _debug:
            print ("_link_prevnode():\n\tfrom=%s\n\tto=%s" %
                   (node, prevnode))

        if prevnode is None:
            return

        node.rev._add_child(prevnode.rev)

    def _unlink_activeleaf(self):
        activeleaf = self.activeleaf
        if activeleaf is None:
            return

        link = activeleaf.link
        if link is None:
            return
        link.rev._delchild(activeleaf.rev)
        activeleaf.link = None

    def _split_fullstr(self):
        # Split the edge leading to the full-string node into a
        # collapsed path.  The new node in the collapsed path becomes
        # a new leaf in the dual tree.
        fullstr = self.fullstr
        if not fullstr:
            return
        parent = fullstr.parent.rev
        if parent.zpath:
            parent = parent.zpath.tail.rev
        wlen = parent._edgelen(fullstr) - 1
        parent._split(fullstr, wlen, True)

    def _update_activeleaf(self, nleaf, term=False):
        """Make activeleaf's link point to newleaf (the most recently
        added leaf).
        """
        activeleaf = self.activeleaf
        if activeleaf is not None:
            nleaf.rev._add_child(activeleaf.rev)

        if not term:
            self.activeleaf = nleaf

    def _term_activeleaf(self):
        # Terminate the active leaf by making a suffix link from it to
        # the explicit node of the new active suffix.
        if not self.isexplicit():
            return
        self._update_activeleaf(self.node, term=True)

    def _update_explicit(self):
        """Update an explicit node reference.

        Return True if finished.
        """
        i = self.stop
        if _debug:
            print ("_update_explicit(%s):\n\t%s" % (i, repr(self)))
        node = self.node
        # Does this ever trigger?  Yes: a split of an edge can be
        # followed by a call to _update_explicit().
        self._link_prevnode(node)
        # Reset prevnode so that the suffix link chain does not
        # collapse, i.e., avoid making improper suffix links directly
        # to the root.
        self.prevnode = None
        if node.haschild(self.tstr[i]):
            # New active suffix is already in the tree.  Walk one
            # character down the existing edge and canonize.
            self.stop += self.step
            self.pathlen += 1
            self.canonize()
            return True
        else:
            # Add leaves, walking suffix links until reaching the root
            # or reaching an explicit node that already contains an
            # edge starting with the character being added.
            if self.activeleaf is None:
                # First leaf represents the whole string, so its label
                # is open at both ends.
                nleaf = node._add_leaf(self.tstr, None)
                self.fullstr = nleaf
            else:
                nleaf = node._add_leaf(self.tstr, i - self.pathlen * self.step)

            self._update_activeleaf(nleaf)
            self.followlink()
            if node.loc is None:
                self.start = self.stop = i + self.step
                return True
            else:
                return False

    def _update_implicit(self):
        """Update an implicit node reference.

        Return True if finished.
        """
        i, node, tstr = self.stop, self.node, self.tstr
        if _debug:
            print ("_update_implicit(%s):\n"
                   "_update_implicit: self=%s" % (i, repr(self)))
        c = tstr[self.start]
        child = node.findchild(c)
        edge = node._edgeword(child)
        wlen = len(self)
        eright = edge.adjstart + wlen * self.step
        if tstr[i] == edge.tstr[eright]:
            # New active suffix is already in the tree.
            self._link_prevnode(node)
            self.stop += self.step
            self.pathlen += 1
            self.canonize()
            return True
        else:
            # Split edges, walking suffix links until reaching the
            # root or an explicit node that already represents the
            # suffix being inserted.  nn is the new node that is
            # created when the edge is split.
            start = i - self.pathlen * self.step
            nn, nleaf = node._split_add_leaf(child, wlen, tstr, start)
            self._link_prevnode(nn)
            self.prevnode = nn
            self._update_activeleaf(nleaf)
            self.followlink()
            return False


class GSTree(SNode, DictMixin):
    """Generalized suffix tree.

    GSTree() -> new empty generalized suffix tree
    GSTree(string, **D) -> new tree from string, using key=0
    GSTree(dict, **D) -> new tree from dict

    In any case, process keyword arguments D as an additional dict of
    strings to add, e.g.,

        t = GSTree(a='abracadabra', b='banana')

    is equivalent to

        t = GSTree({'a': 'abracadabra', 'b': 'banana'})

    >>> s = 'amanaplanacanalpanamazz'
    >>> t = GSTree(a=s, b=s[::-1])
    >>> t.printpairs(minlen=4)
    a[0:21] b[2:23] amanaplanacanalpanama

    GSTree objects behave somewhat like dictionaries whose values are
    the complete strings contained in the tree.  Complete strings may
    be added using mapping syntax, e.g.

    >>> t = GSTree(a='eh')
    >>> t['foo'] = 'foo1'
    >>> t['bar'] = 'blah'

    but retrieving the strings results in a wrapped object that may
    implement mutation in the future.

    >>> str(t['a'])
    'eh'
    >>> str(t['bar'])
    'blah'
    """
    def __init__(self, _arg=None, _nobuild=False, **kwargs):
        self._strings = {}
        self._nextkey = 0
        self._nobuild = _nobuild
        super(GSTree, self).__init__()
        if _arg is not None:
            if hasattr(_arg, 'keys'):
                self.addstrings(_arg)
            else:
                self.addstring(_arg)
        self.addstrings(**kwargs)

    def __repr__(self):
        d = dict((k, v.s) for (k, v) in self._strings.iteritems())
        return "%s(%s)" % (self.__class__.__name__, d)

    def __setitem__(self, key, s):
        self._addstring(key, s)

    def __getitem__(self, key):
        return self._strings[key]

    def keys(self):
        return self._strings.keys()

    # Map key to terminator.  Override this in a subclass if the
    # alphabet can have variable-length characters or other weirdness.
    def _key2term(self, key):
        return '$' + str(key)

    # Map terminator to key.  Return None if not a terminator.
    # Override this in a subclass if the alphabet can have
    # variable-length characters or other weirdness.
    def _term2key(self, term):
        if len(term) > 1 and term[0] == '$':
            return term[1:]
        else:
            return None

    def _key2init(self, key):
        return '^' + str(key)

    def _init2key(self, init):
        if len(init) > 1 and init[0] == '^':
            return init[1:]
        else:
            return None

    def chr(self, key, i):
        """The character at position i of the string with the specified key."""
        return self._strings[key][i]

    def _addstring(self, key, s):
        if key in self._strings:
            raise ValueError("key %s already exists in tree" % key)
        if hasattr(s, '__setitem__'):
            raise TypeError('string must be immutable')
        tstr = TreeStr(self, key, s)
        self._strings[key] = tstr

        if self._nobuild:
            return

        active = ActiveRef(self, tstr, 0, 0, 1)
        for i in xrange(len(s) + 1):
            active.update(i)

    def addstring(self, *args):
        """
        Add a string to the tree.

        t.addstring(string) -> add string to tree using next available key
        t.addstring(key, string) -> add string to tree using specified key
        """
        if len(args) == 1:
            s = args[0]
            key = self._nextkey
            self._nextkey += 1
        elif len(args) == 2:
            key, s = args
        else:
            raise TypeError('%s.addstring: invalid number of arguments' %
                            (self.__class__.__name__,))
        self._addstring(key, s)

    def addstrings(self, *args, **kwargs):
        """
        Add multiple strings to the tree.

        t.addstrings(dict, **D) -> update tree from dict
        t.addstrings(str1, str2, ..., **D) -> add strings to tree

        In any case, process keyword arguments **D as an additional
        dict of strings to add.
        """
        if len(args) == 1:
            arg = args[0]
            if hasattr(arg, 'keys'):
                for k in arg:
                    self.addstring(k, arg[k])
            else:
                self.addstring(arg)
        else:
            for s in args:
                self.addstring(s)
        for k, v in kwargs.iteritems():
            self.addstring(k, v)

    def _sortkey(self, c):
        """Sort key for string characters, accounting for terminators.
        Terminators compare less than all string characters.
        """
        tkey = self._term2key(c)
        if tkey is None:
            return (1, c)
        else:
            return (0, tkey)

    def _paircoords(self, pair):
        ((lkey, lleft), (rkey, rleft), wordlen) = pair
        lright = lleft + wordlen
        rright = rleft + wordlen
        return ("%s[%s:%s] %s[%s:%s]" %
                (lkey, lleft, lright, rkey, rleft, rright))

    def _pairword(self, pair):
        ((key, left), (k2, l2), wordlen) = pair
        right = left + wordlen
        return self._strings[key][left:right]

    def printpairs(self, pairs=None, **kwargs):
        """Print a list of pairs (as from maxpairs() method) in a more
        readable format.  Output looks like:

            key1[start1:stop1] key2[start2:stop2] word

        If no pair list is given, run self.itermaxpairs() and format its
        output.

        >>> t = GSTree('mississippi')
        >>> t.printpairs(sorted(t.itermaxpairs()))
        0[1:5] 0[4:8] issi
        0[1:2] 0[7:8] i
        0[1:2] 0[10:11] i
        0[2:3] 0[3:4] s
        0[2:3] 0[6:7] s
        0[3:4] 0[5:6] s
        0[4:5] 0[10:11] i
        0[5:6] 0[6:7] s
        0[7:8] 0[10:11] i
        0[8:9] 0[9:10] p
        """
        if pairs is None:
            pairs = self.itermaxpairs(**kwargs)
        for p in pairs:
            print "%s %s" % (self._paircoords(p), self._pairword(p))


class IterNodes(object):
    """Iterator that recursively walks a subtree.

    IterNodes(node) -> iterator over the subtree rooted at node, not an actual
        instance object.  Effectively the same as IterNodes()(node).

    IterNodes() -> instance

    Instances are callables that take a starting node as parameter.
    When called, an instance walks the subtree rooted at the specified
    node.  Subclasses can override the callable behavior, e.g., making
    it take more parameters.
    """
    def __new__(cls, *args, **kwargs):
        inst = super(IterNodes, cls).__new__(cls)
        if args or kwargs:
            return inst(*args, **kwargs)
        else:
            return inst

    def __call__(self, node):
        for x in self._step([node]):
            yield x

    def _children(self, nodes):
        return nodes[-1].iterchildren()

    def _sortchildren(self, nodes):
        return nodes[-1].sortedchildren()

    def _noproc(self, nodes):
        return ()

    def _yesproc(self, nodes):
        return (nodes,)

    def _step(self, nodes):
        for x in self.preproc(nodes):
            yield x
        for child in self.childproc(nodes):
            for x in self._step(nodes + [child]):
                yield x
        for x in self.postproc(nodes):
            yield x

    childproc = _children
    preproc = _yesproc
    postproc = _noproc


class IterNodeStrs(IterNodes):
    """Iterator that recursively produces string representations of
    nodes, including the paths leading to them.
    """
    def __call__(self, node, hexids=True):
        self.hexids = hexids
        return super(IterNodeStrs, self).__call__(node)

    def _pathstr(self, nodes):
        g = izip(nodes, islice(nodes, 1, None))
        s = ' '.join(parent._edgeword(child).strword() for parent, child in g)
        try:
            s += " (len=%d)" % len(nodes[-1].loc)
        except TypeError:
            pass
        if self.hexids:
            s += "%s%s" % (' ' if s else '', hex(id(nodes[-1])))
        return s

    def _nodestr(self, nodes):
        if not self.leavesonly or not nodes[-1].child:
            yield self._pathstr(nodes)

    leavesonly = False
    childproc = IterNodes._sortchildren

    preproc = _nodestr


class IterLeafStrs(IterNodeStrs):
    """Like IterNodeStrs, but restricted to leaf nodes."""
    leavesonly = True


class IterWhere(IterNodes):
    """Iterator looks for a word in a tree (or subtree).

    Instances return all locations of the specified word, as a list of
    (key, location) tuples.
    """
    def __call__(self, node, word):
        self.word = word
        return super(IterWhere, self).__call__(node)

    def _where_filt(self, nodes):
        word = self.word
        node = nodes[-1]
        # print "_where_filt: word=%s" % (word,)
        # _printnode(walker, node, parents)
        if len(word) == 0:
            return node.iterchildren()
        child = node.findchild(word[0])
        if child is None:
            return ()
        elen = node._edgelen(child)
        wlen = len(word)
        eword = node._edgeword(child).word()
        if word[0:elen] == eword:
            self.word = word[elen:]
        elif eword[0:wlen] == word:
            self.word = word[wlen:]
        else:
            return ()
        return (child,)

    def _where_step(self, nodes):
        if len(nodes) <= 1:
            return
        node = nodes[-1]
        if len(self.word) != 0:
            return
        if node.isleaf():
            loc = node.loc
            yield (loc.tstr.key, loc.adjstart)

    childproc = _where_filt
    preproc = _where_step


class IterMaxPairs(IterNodes):
    """IterNodes subclass that iterates over maximal pairs."""
    def __call__(self, node, minlen=0, allpairs=False):
        self.minlen = minlen
        self.allpairs = allpairs
        return super(IterMaxPairs, self).__call__(node)

    def _maxpair_step(self, nodes):
        if len(nodes) <= 1:
            return
        node = nodes[-1]
        if len(node.loc) < self.minlen:
            return
        if node.isleaf():
            self._maxpair_leaf(node)
            return
        else:
            edgelen = nodes[-2]._edgelen(nodes[-1])
            for pair in self._maxpair_branch(node, edgelen):
                yield pair

    def _maxpair_leaf(self, node):
        # A leaf node denotes a suffix.  Determine the character (the
        # lchr) to the left of (immediately preceding) the suffix,
        # along with its position.
        loc = node.loc
        pos = loc.adjstart
        lpos = pos - loc.step
        if pos <= 0:
            # Unique per key, using initiator, otherwise pairs where
            # each occurrence is at the beginning of a different
            # string will not be found.
            lchr = loc.tstr.init
        else:
            lchr = loc.tstr[lpos]
        node.lchr = lchr
        # Initialize the lchrlist for this lchr.  A leaf node has
        # exactly one lchr and lchr position.
        node.lchrlists[lchr] = _list(((loc.tstr.key, pos),))

    def _maxpair_branch(self, node, edgelen):
        pathlen = len(node.loc)
        for lchild, rchild in combinations(node.iterchildren(), 2):
            iterlc = lchild.lchrlists.iteritems()
            iterrc = rchild.lchrlists.iteritems()
            for poslpair in product(iterlc, iterrc):
                for pair in self._maxpair_child(poslpair, pathlen, edgelen):
                    yield pair

        for child in node.iterchildren():
            for lchr, posl in child.lchrlists.iteritems():
                # Consolidate the lchrlists of the children into this
                # node.  The copying might be slower than consing
                # linked lists, if the lists are large.
                if lchr not in node.lchrlists:
                    node.lchrlists[lchr] = _list()
                node.lchrlists[lchr].extend(posl)

    def _maxpair_child(self, poslpair, pathlen, edgelen):
        minlen = self.minlen
        allpairs = self.allpairs
        (llchr, lposl), (rlchr, rposl) = poslpair
        if llchr == rlchr and not allpairs:
            return
        for lpos, rpos in product(lposl, rposl):
            left, right = sorted((lpos, rpos))
            if not allpairs:
                yield (left, right, pathlen)
            else:
                for i in xrange(minlen, pathlen + 1):
                    yield (left, right, i)

    def _clear_lchrs(self, nodes):
        nodes[-1].lchrlists = {}
        return ()

    preproc = _clear_lchrs
    postproc = _maxpair_step


def _randstr():
    import string
    import random
    return ''.join(random.choice(string.letters) for _ in xrange(5000))

def _randdna():
    import random
    return ''.join(random.choice('ATGC') for _ in xrange(5000))

def _timetest():
    import timeit

    print 'Timing random letters...'
    out = timeit.repeat('GSTree(s)', """\
from %s import _randstr, GSTree
s = _randstr()""" % __name__,
                        repeat=3, number=10)
    print out

    print 'Timing random DNA sequences...'
    out = timeit.repeat('GSTree(s)', """\
from %s import _randdna, GSTree
s = _randdna()""" % __name__,
                        repeat=3, number=10)
    print out

if __name__ == '__main__':
    _timetest()
