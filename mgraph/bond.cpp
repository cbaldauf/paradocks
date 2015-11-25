#include "bond.hpp"

using namespace boost::assign;

namespace mgraph 
{

const string BT::_convertbt2str[6] = {"s", "re", "ar", "d", "t", "UNK"};

const map<string, BT::_BT> BT::_convertstr2bt = map_list_of
("1",  BT::s)  ("s", BT::s)
("re", BT::re)
("ar", BT::ar)
("2",  BT::d)  ("d", BT::d)
("3",  BT::t)  ("t", BT::t)
("UNK", BT::UNK);

}
